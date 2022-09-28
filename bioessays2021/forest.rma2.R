#### HIDDEN FUNCTIONS ####

############################################################################

### function to set default 'btt' value(s) or check specified 'btt' values

.set.btt <- function(btt, p, int.incl, X) {
  
  mstyle <- .get.mstyle("crayon" %in% .packages())
  
  if (missing(btt) || is.null(btt)) {
    
    if (p > 1) {                        ### if the model matrix has more than one column
      if (int.incl) {
        btt <- seq.int(from=2, to=p)     ### and the model has an intercept term, test all coefficients except the intercept
      } else {
        btt <- seq_len(p)                ### and the model does not have an intercept term, test all coefficients
      }
    } else {
      btt <- 1                         ### if the model matrix has a single column, test that single coefficient
    }
    
  } else {
    
    if (is.character(btt)) {
      
      btt <- grep(btt, colnames(X))
      
      if (length(btt) == 0L)
        stop(mstyle$stop("Cannot identify coefficient(s) corresponding to the specified 'btt' string."))
      
    } else {
      
      ### round, take unique values, and sort
      btt <- sort(unique(round(btt)))
      
      ### check for mix of positive and negative values
      if (any(btt < 0) && any(btt > 0))
        stop(mstyle$stop("Cannot mix positive and negative 'btt' values."))
      
      ### keep/remove from 1:p vector as specified
      btt <- seq_len(p)[btt]
      
      ### (1:5)[5:6] yields c(5, NA) so remove NAs if this happens
      btt <- btt[!is.na(btt)]
      
      ### make sure that at least one valid value is left
      if (length(btt) == 0L)
        stop(mstyle$stop("Non-existent coefficients specified via 'btt'."))
      
    }
    
  }
  
  return(btt)
  
}

### function to format 'btt' values for printing

.format.btt <- function(btt) {
  
  sav <- c()
  
  if (length(btt) > 1L) {
    
    while (length(btt) > 0L) {
      
      x <- rle(diff(btt))
      
      if (x$values[1] == 1 && length(x$values) != 0L) {
        sav <- c(sav, c(btt[1], ":", btt[x$lengths[1] + 1]))
        btt <- btt[-c(1:(x$lengths[1] + 1))]
        sav <- c(sav, ", ")
      } else {
        sav <- c(sav, btt[1], ",")
        btt <- btt[-1]
      }
      
    }
    
    sav <- paste0(sav[-length(sav)], collapse="")
    
  } else {
    
    sav <- paste0(btt)
    
  }
  
  return(sav)
  
}

############################################################################

### pairwise sorting of the elements of two vectors

.psort <- function(x,y) {
  
  ### t(apply(xy, 1, sort)) would be okay, but problematic if there are NAs;
  ### either they are removed completely (na.last=NA) or they are always put
  ### first/last (na.last=FALSE/TRUE); but we just want to leave the NAs in
  ### their position!
  
  if (is.null(x) || length(x) == 0L) ### need to catch this
    return(NULL)
  
  if (missing(y)) {
    if (is.matrix(x)) {
      xy <- x
    } else {
      xy <- rbind(x) ### in case x is just a vector
    }
  } else {
    xy <- cbind(x,y)
  }
  
  n <- nrow(xy)
  
  for (i in seq_len(n)) {
    if (anyNA(xy[i,]))
      next
    xy[i,] <- sort(xy[i,])
  }
  
  colnames(xy) <- NULL
  
  return(xy)
  
}

############################################################################

### function to obtain the trace of a matrix

.tr <- function(X)
  return(sum(diag(X)))

### function to check if a matrix is square

.is.square <- function(X)
  NROW(X) == NCOL(X)

### use NROW/NCOL to better deal with scalars; compare:
### (V <- list(matrix(1, nrow=2, ncol=2), 3, c(1,4), cbind(c(2,1)))); sapply(V, function(x) nrow(x) == ncol(x)); sapply(V, function(x) NROW(x) == NCOL(x))

### function to test whether a vector is all equal to 1s (e.g., to find intercept(s) in a model matrix)

.is.intercept <- function(x, eps=1e-08)
  return(all(abs(x - 1) < eps))

### function to test whether a vector is a dummy variable (i.e., consists of only 0s and 1s)

.is.dummy <- function(x, eps=1e-08)
  return(all(abs(x) < eps | abs(x - 1) < eps))
#return(all(sapply(x, identical, 0) | sapply(x, identical, 1)))

### function to test whether something is a vector (in the sense of being atomic, not a matrix, and not NULL)

.is.vector <- function(x)
  is.atomic(x) && !is.matrix(x) && !is.null(x)

############################################################################

### function to format p-values
### if showeq=FALSE, c(.001, .00001) becomes c("0.0010", "<.0001")
### if showeq=TRUE,  c(.001, .00001) becomes c("=0.0010", "<.0001")
### if add0=FALSE, "<.0001"; if add0=TRUE, "<0.0001"

.pval <- function(p, digits=4, showeq=FALSE, sep="", add0=FALSE) {
  
  digits <- max(digits, 1)
  cutoff  <- paste(c(".", rep(0,digits-1),1), collapse="")
  ncutoff <- as.numeric(cutoff)
  
  ifelse(is.na(p), paste0(ifelse(showeq, "=", ""), sep, NA),
         ifelse(p >= ncutoff, paste0(ifelse(showeq, "=", ""), sep, formatC(p, digits=digits, format="f")),
                paste0("<", sep, ifelse(add0, "0", ""), cutoff)))
  
}

### function to format/round values in general

.fcf <- function(x, digits) {
  
  if (all(is.na(x))) { # since formatC(NA, format="f", digits=2) fails
    x
  } else {
    formatC(x, format="f", digits=digits)
  }
  
}

############################################################################

### function to print a named (character) vector right aligned with
### a gap of two spaces between adjacent values and no padding

.print.vector <- function(x) {
  
  if (is.null(names(x)))
    names(x) <- seq_along(x)
  
  len.n   <- nchar(names(x))
  len.x   <- nchar(x)
  len.max <- pmax(len.n, len.x)
  format  <- sapply(len.max, function(x) paste("%", x, "s", sep=""))
  
  row.n <- paste(sprintf(format, names(x)), collapse="  ")
  row.x <- paste(sprintf(format, x), collapse="  ")
  
  cat(row.n, "\n", row.x, "\n", sep="")
  
}

############################################################################

### function like make.unique(), but starts at .1 for the first instance
### of a repeated element

.make.unique <- function(x) {
  
  x <- as.character(x)
  ux <- unique(x)
  
  for (i in seq_along(ux)) {
    xiTF <- x == ux[i]
    xi <- x[xiTF]
    if (length(xi) == 1L)
      next
    x[xiTF] <- paste(xi, seq_along(xi), sep=".")
  }
  
  return(x)
  
}

############################################################################

### function to check if extra/superfluous arguments are specified via ...

.chkdots <- function(ddd, okargs) {
  
  mstyle <- .get.mstyle("crayon" %in% .packages())
  
  for (i in seq_along(okargs))
    ddd[okargs[i]] <- NULL
  
  if (length(ddd) > 0L)
    warning(mstyle$warning(paste0("Extra argument", ifelse(length(ddd) > 1L, "s ", " "), "(", paste0("'", names(ddd), "'", collapse=", "), ") disregarded.")), call.=FALSE)
  
}

############################################################################

### set axis label (for forest, funnel, and labbe functions)

.setlab <- function(measure, transf.char, atransf.char, gentype, short=FALSE) {
  
  if (gentype == 1)
    lab <- "Observed Outcome"
  if (gentype == 2)
    lab <- "Overall Estimate" # for forest.cumul.rma() function
  if (gentype == 3)
    lab <- "Estimate"         # for header
  
  #########################################################################
  
  if (!is.null(measure)) {
    
    ######################################################################
    if (is.element(measure, c("RR","MPRR"))) {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[RR]", "Log Risk Ratio")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Risk Ratio")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Risk Ratio", "Risk Ratio (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Risk Ratio", "Risk Ratio")
      }
    }
    if (is.element(measure, c("OR","PETO","D2OR","D2ORN","D2ORL","MPOR","MPORC","MPPETO"))) {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[OR]", "Log Odds Ratio")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Odds Ratio")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Odds Ratio", "Odds Ratio (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Odds Ratio", "Odds Ratio")
      }
    }
    if (is.element(measure, c("RD","MPRD"))) {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Risk Difference", "Risk Difference")
      } else {
        lab <- ifelse(short, lab, "Transformed Risk Difference")
      }
    }
    if (measure == "AS") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Arcsine RD", "Arcsine Transformed Risk Difference")
      } else {
        lab <- ifelse(short, lab, "Transformed Arcsine Transformed Risk Difference")
      }
    }
    if (measure == "PHI") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Phi", "Phi Coefficient")
      } else {
        lab <- ifelse(short, lab, "Transformed Phi Coefficient")
      }
    }
    if (measure == "YUQ") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Yule's Q", "Yule's Q")
      } else {
        lab <- ifelse(short, lab, "Transformed Yule's Q")
      }
    }
    if (measure == "YUY") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Yule's Y", "Yule's Y")
      } else {
        lab <- ifelse(short, lab, "Transformed Yule's Y")
      }
    }
    ######################################################################
    if (measure == "IRR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[IRR]", "Log Incidence Rate Ratio")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Incidence Rate Ratio")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Rate Ratio", "Incidence Rate Ratio (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Rate Ratio", "Incidence Rate Ratio")
      }
    }
    if (measure == "IRD") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "IRD", "Incidence Rate Difference")
      } else {
        lab <- ifelse(short, lab, "Transformed Incidence Rate Difference")
      }
    }
    if (measure == "IRSD") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "IRSD", "Square Root Transformed Incidence Rate Difference")
      } else {
        lab <- ifelse(short, lab, "Transformed Square Root Transformed Incidence Rate Difference")
      }
    }
    ######################################################################
    if (measure == "MD") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "MD", "Mean Difference")
      } else {
        lab <- ifelse(short, lab, "Transformed Mean Difference")
      }
    }
    if (is.element(measure, c("SMD","SMDH","PBIT","OR2D","OR2DN","OR2DL"))) {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "SMD", "Standardized Mean Difference")
      } else {
        lab <- ifelse(short, lab, "Transformed Standardized Mean Difference")
      }
    }
    if (measure == "ROM") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[RoM]", "Log Ratio of Means")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Ratio of Means")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Ratio of Means", "Ratio of Means (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Ratio of Means", "Ratio of Means")
      }
    }
    if (measure == "RPB") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Correlation", "Point-Biserial Correlation")
      } else {
        lab <- ifelse(short, lab, "Transformed Point-Biserial Correlation")
      }
    }
    if (measure == "CVR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[CVR]", "Log Coefficient of Variation Ratio")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation Ratio")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio")
      }
    }
    if (measure == "VR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[VR]", "Log Variability Ratio")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Variability Ratio")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "VR", "Variability Ratio (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "VR", "Variability Ratio")
      }
    }
    ######################################################################
    if (is.element(measure, c("COR","UCOR","RTET","RBIS"))) {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Correlation", "Correlation Coefficient")
      } else {
        lab <- ifelse(short, lab, "Transformed Correlation Coefficient")
      }
    }
    if (measure == "ZCOR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, expression('Fisher\'s ' * z[r]), "Fisher's z Transformed Correlation Coefficient")
      } else {
        lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Correlation Coefficient")
        if (atransf.char == "transf.ztor" || atransf.char == "transf.ztor.int")
          lab <- ifelse(short, "Correlation", "Correlation Coefficient")
        if (transf.char == "transf.ztor" || transf.char == "transf.ztor.int")
          lab <- ifelse(short, "Correlation", "Correlation Coefficient")
      }
    }
    ######################################################################
    if (measure == "PCOR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
      } else {
        lab <- ifelse(short, lab, "Transformed Partial Correlation Coefficient")
      }
    }
    if (measure == "ZPCOR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, expression('Fisher\'s ' * z[r]), "Fisher's z Transformed Partial Correlation Coefficient")
      } else {
        lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Partial Correlation Coefficient")
        if (atransf.char == "transf.ztor" || atransf.char == "transf.ztor.int")
          lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
        if (transf.char == "transf.ztor" || transf.char == "transf.ztor.int")
          lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
      }
    }
    if (measure == "SPCOR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Correlation", "Semi-Partial Correlation Coefficient")
      } else {
        lab <- ifelse(short, lab, "Transformed Semi-Partial Correlation Coefficient")
      }
    }
    ######################################################################
    if (measure == "PR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Proportion", "Proportion")
      } else {
        lab <- ifelse(short, lab, "Transformed Proportion")
      }
    }
    if (measure == "PLN") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[Pr]", "Log Proportion")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Proportion")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Proportion", "Proportion (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Proportion", "Proportion")
      }
    }
    if (measure == "PLO") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[Odds]", "Log Odds")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Odds")
        if (atransf.char == "transf.ilogit" || atransf.char == "transf.ilogit.int" || atransf.char == "plogis")
          lab <- ifelse(short, "Proportion", "Proportion (logit scale)")
        if (transf.char == "transf.ilogit" || transf.char == "transf.ilogit.int" || transf.char == "plogis")
          lab <- ifelse(short, "Proportion", "Proportion")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Odds", "Odds (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Odds", "Odds")
      }
    }
    if (measure == "PAS") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, expression(arcsin(sqrt(p))), "Arcsine Transformed Proportion")
      } else {
        lab <- ifelse(short, lab, "Transformed Arcsine Transformed Proportion")
        if (atransf.char == "transf.iarcsin" || atransf.char == "transf.iarcsin.int")
          lab <- ifelse(short, "Proportion", "Proportion (arcsine scale)")
        if (transf.char == "transf.iarcsin" || transf.char == "transf.iarcsin.int")
          lab <- ifelse(short, "Proportion", "Proportion")
      }
    }
    if (measure == "PFT") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "PFT", "Double Arcsine Transformed Proportion")
      } else {
        lab <- ifelse(short, lab, "Transformed Double Arcsine Transformed Proportion")
        if (atransf.char == "transf.ipft.hm")
          lab <- ifelse(short, "Proportion", "Proportion")
        if (transf.char == "transf.ipft.hm")
          lab <- ifelse(short, "Proportion", "Proportion")
      }
    }
    ######################################################################
    if (measure == "IR") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Rate", "Incidence Rate")
      } else {
        lab <- ifelse(short, lab, "Transformed Incidence Rate")
      }
    }
    if (measure == "IRLN") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[IR]", "Log Incidence Rate")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Incidence Rate")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Rate", "Incidence Rate (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Rate", "Incidence Rate")
      }
    }
    if (measure == "IRS") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Sqrt[IR]", "Square Root Transformed Incidence Rate")
      } else {
        lab <- ifelse(short, lab, "Transformed Square Root Transformed Incidence Rate")
        if (atransf.char == "transf.isqrt" || atransf.char == "transf.isqrt.int")
          lab <- ifelse(short, "Rate", "Incidence Rate (square root scale)")
        if (transf.char == "transf.isqrt" || transf.char == "transf.isqrt.int")
          lab <- ifelse(short, "Rate", "Incidence Rate")
      }
    }
    if (measure == "IRFT") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "IRFT", "Freeman-Tukey Transformed Incidence Rate")
      } else {
        lab <- ifelse(short, lab, "Transformed Freeman-Tukey Transformed Incidence Rate")
      }
    }
    ######################################################################
    if (measure == "MN") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Mean", "Mean")
      } else {
        lab <- ifelse(short, lab, "Transformed Mean")
      }
    }
    if (measure == "MNLN") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[Mean]", "Log Mean")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Mean")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Mean", "Mean (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Mean", "Mean")
      }
    }
    if (measure == "CVLN") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[CV]", "Log Coefficient of Variation")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "CV", "Coefficient of Variation (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "CV", "Coefficient of Variation")
      }
    }
    if (measure == "SDLN") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[SD]", "Log Standard Deviation")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Standard Deviation")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "SD", "Standard Deviation (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "SD", "Standard Deviation")
      }
    }
    ######################################################################
    if (measure == "MC") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Mean Change", "Mean Change")
      } else {
        lab <- ifelse(short, lab, "Transformed Mean Change")
      }
    }
    if (is.element(measure, c("SMCC","SMCR","SMCRH"))) {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "SMC", "Standardized Mean Change")
      } else {
        lab <- ifelse(short, lab, "Transformed Standardized Mean Change")
      }
    }
    if (measure == "ROMC") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[RoM]", "Log Ratio of Means")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Ratio of Means")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "Ratio of Means", "Ratio of Means (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "Ratio of Means", "Ratio of Means")
      }
    }
    if (measure == "CVRC") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[CVR]", "Log Coefficient of Variation Ratio")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation Ratio")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio")
      }
    }
    if (measure == "VRC") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Log[VR]", "Log Variability Ratio")
      } else {
        lab <- ifelse(short, lab, "Transformed Log Variability Ratio")
        if (atransf.char == "exp" || atransf.char == "transf.exp.int")
          lab <- ifelse(short, "VR", "Variability Ratio (log scale)")
        if (transf.char == "exp" || transf.char == "transf.exp.int")
          lab <- ifelse(short, "VR", "Variability Ratio")
      }
    }
    ######################################################################
    if (measure == "ARAW") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, "Alpha", "Cronbach's alpha")
      } else {
        lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
      }
    }
    if (measure == "AHW") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, expression('Alpha'[HW]), "Transformed Cronbach's alpha")
      } else {
        lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
        if (atransf.char == "transf.iahw")
          lab <- ifelse(short, "Alpha", "Cronbach's alpha")
        if (transf.char == "transf.iahw")
          lab <- ifelse(short, "Alpha", "Cronbach's alpha")
      }
    }
    if (measure == "ABT") {
      if (transf.char == "FALSE" && atransf.char == "FALSE") {
        lab <- ifelse(short, expression('Alpha'[B]), "Transformed Cronbach's alpha")
      } else {
        lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
        if (atransf.char == "transf.iabt")
          lab <- ifelse(short, "Alpha", "Cronbach's alpha")
        if (transf.char == "transf.iabt")
          lab <- ifelse(short, "Alpha", "Cronbach's alpha")
      }
    }
    ######################################################################
    
  }
  
  return(lab)
  
}

############################################################################

### stuff related to colored/styled output

.get.mstyle <- function(withcrayon) {
  
  if (withcrayon) {
    
    if (exists(".mstyle")) {
      .mstyle <- get(".mstyle")
      if (!is.list(.mstyle))
        .mstyle <- list(.mstyle)
    } else {
      .mstyle <- list()
    }
    
    if (is.null(.mstyle$section)) {
      section <- crayon::bold
    } else {
      section <- .mstyle$section
    }
    if (is.null(.mstyle$header)) {
      header <- crayon::underline
    } else {
      header <- .mstyle$header
    }
    if (is.null(.mstyle$body)) {
      body <- crayon::reset
    } else {
      body <- .mstyle$body
    }
    if (is.null(.mstyle$text)) {
      text <- crayon::reset
    } else {
      text <- .mstyle$text
    }
    if (is.null(.mstyle$result)) {
      result <- crayon::reset
    } else {
      result <- .mstyle$result
    }
    if (is.null(.mstyle$stop)) {
      stop <- crayon::combine_styles(crayon::red, crayon::bold)
    } else {
      stop <- .mstyle$stop
    }
    if (is.null(.mstyle$warning)) {
      warning <- crayon::yellow
    } else {
      warning <- .mstyle$warning
    }
    if (is.null(.mstyle$message)) {
      message <- crayon::green
    } else {
      message <- .mstyle$message
    }
    if (is.null(.mstyle$verbose)) {
      verbose <- crayon::cyan
    } else {
      verbose <- .mstyle$verbose
    }
    if (is.null(.mstyle$legend)) {
      legend <- crayon::silver
    } else {
      legend <- .mstyle$legend
    }
    
  } else {
    
    tmp <- function(...) paste0(...)
    section <- tmp
    header  <- tmp
    body    <- tmp
    text    <- tmp
    result  <- tmp
    stop    <- tmp
    warning <- tmp
    message <- tmp
    verbose <- tmp
    legend  <- tmp
    
  }
  
  return(list(section=section, header=header, body=body, text=text, result=result, stop=stop, warning=warning, message=message, verbose=verbose, legend=legend))
  
}

.print.output <- function(x, mstyle) {
  
  if (missing(mstyle)) {
    for (i in seq_along(x)) {
      cat(x[i], "\n")
    }
  } else {
    for (i in seq_along(x)) {
      cat(mstyle(x[i]), "\n")
    }
  }
  
}

.print.table <- function(x, mstyle) {
  
  is.header <- !grepl(" [-0-9]", x)
  
  for (i in seq_along(x)) {
    if (is.header[i]) {
      x[i] <- trimws(x[i], which="right")
      x[i] <- mstyle$header(x[i])
    } else {
      x[i] <- mstyle$body(x[i])
    }
    cat(x[i], "\n")
  }
  
}

############################################################################

.set.digits <- function(digits, dmiss) {
  
  res <- c(est=4, se=4, test=4, pval=4, ci=4, var=4, sevar=4, fit=4, het=4)
  
  if (exists(".digits")) {
    .digits <- get(".digits")
    if (is.null(names(.digits)) && length(.digits) == 1L) {
      # if .digits is a single unnamed scalar, set all digit values to that value
      res <- c(est=.digits, se=.digits, test=.digits, pval=.digits, ci=.digits, var=.digits, sevar=.digits, fit=.digits, het=.digits)
    } else if (any(names(.digits) != "") && any(names(.digits) == "")) {
      # if .digits has (at least) one unnamed element, use it to set all unnamed elements to that digits value
      pos <- pmatch(names(.digits), names(res))
      res[c(na.omit(pos))] <- .digits[!is.na(pos)]
      otherval <- .digits[names(.digits) == ""][1]
      res[(1:9)[-c(na.omit(pos))]] <- otherval
    } else {
      pos <- pmatch(names(.digits), names(res))
      res[c(na.omit(pos))] <- .digits[!is.na(pos)]
    }
  }
  
  if (!dmiss) {
    if (is.null(names(digits))) {
      res <- c(est=digits[[1]], se=digits[[1]], test=digits[[1]], pval=digits[[1]], ci=digits[[1]], var=digits[[1]], sevar=digits[[1]], fit=digits[[1]], het=digits[[1]])
    } else {
      pos <- pmatch(names(digits), names(res))
      res[c(na.omit(pos))] <- digits[!is.na(pos)]
    }
  }
  
  res
  
}

.get.digits <- function(digits, xdigits, dmiss) {
  
  res <- xdigits
  
  if (exists(".digits")) {
    .digits <- get(".digits")
    pos <- pmatch(names(.digits), names(res))
    res[c(na.omit(pos))] <- .digits[!is.na(pos)]
  }
  
  if (!dmiss) {
    if (is.null(names(digits))) {
      res <- c(est=digits[[1]], se=digits[[1]], test=digits[[1]], pval=digits[[1]], ci=digits[[1]], var=digits[[1]], sevar=digits[[1]], fit=digits[[1]], het=digits[[1]])
    } else {
      pos <- pmatch(names(digits), names(res))
      res[c(na.omit(pos))] <- digits[!is.na(pos)]
    }
  }
  
  ### so we can still print objects created with older metafor versions (where xdigit will be just an unnamed scalar)
  if (length(res) == 1L && is.null(names(res)))
    res <- c(est=res[[1]], se=res[[1]], test=res[[1]], pval=res[[1]], ci=res[[1]], var=res[[1]], sevar=res[[1]], fit=res[[1]], het=res[[1]])
  
  res
}

############################################################################

### check if x is logical and TRUE/FALSE (NAs and NULL always evaluate as FALSE)

.isTRUE <- function(x)
  !is.null(x) && is.logical(x) && !is.na(x) && x

.isFALSE <- function(x)
  !is.null(x) && is.logical(x) && !is.na(x) && !x

############################################################################

### to register getfit method for 'rma.uni' and 'rma.mv' objects: eval(metafor:::.glmulti)

.glmulti <- parse(text="

if (!(\"glmulti\" %in% .packages()))
   stop(\"Need to load the 'glmulti' package first to use this code.\")

setOldClass(\"rma.uni\")

setMethod(\"getfit\", \"rma.uni\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

setOldClass(\"rma.mv\")

setMethod(\"getfit\", \"rma.mv\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

setOldClass(\"rma.glmm\")

setMethod(\"getfit\", \"rma.glmm\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

")

### helper functions to make MuMIn work together with metafor

.MuMIn <- parse(text="

makeArgs.rma <- function (obj, termNames, comb, opt, ...) {
   ret <- MuMIn:::makeArgs.default(obj, termNames, comb, opt)
   names(ret)[1L] <- \"mods\"
   ret
}

coefTable.rma <- function (model, ...) {
  MuMIn:::.makeCoefTable(model$b, model$se, coefNames = rownames(model$b))
}

")

### helper functions to make mice work together with metafor

.mice <- parse(text="

glance.rma <- function (x, ...)
   data.frame(df.residual=df.residual(x))

tidy.rma <- function (x, ...) {
   ret <- coef(summary(x))
   colnames(ret)[2] <- \"std.error\"
   ret$term <- rownames(ret)
   return(ret)
}

")

############################################################################

### shorten a string vector so that elements remain distinguishable

.shorten <- function(x, minlen) {
  
  y <- x
  
  x <- c(na.omit(x))
  
  n <- length(unique(x))
  
  maxlen <- max(nchar(unique(x)))
  
  for (l in 1:maxlen) {
    tab <- table(x, substr(x, 1, l))
    if (nrow(tab) == n && ncol(tab) == n && sum(tab[upper.tri(tab)]) == 0 && sum(tab[lower.tri(tab)]) == 0)
      break
  }
  
  if (!missing(minlen) && l < minlen) {
    if (minlen > maxlen)
      minlen <- maxlen
    l <- minlen
  }
  
  return(substr(y, 1, l))
  
}

############################################################################

### simplified version of what mvtnorm::rmvnorm() does

.mvrnorm <- function(n, mu, Sigma) {
  
  p <- nrow(Sigma)
  eS <- eigen(Sigma, symmetric = TRUE)
  eval <- eS$values
  evec <- eS$vectors
  
  Y <- matrix(rnorm(p * n), nrow = n, byrow = TRUE) %*% t(evec %*% (t(evec) * sqrt(pmax(eval, 0))))
  Y <- sweep(Y, 2, mu, "+")
  
  return(Y)
  
}

############################################################################


#### NEW FOREST.RMA(...) ####

forest.rma2 <- function(x, annotate=TRUE, addfit=TRUE, addcred=FALSE, 
                        showweights=FALSE, header=FALSE,
                        xlim, alim, clim, ylim, top=3, at, steps=5, level=x$level, 
                        refline=0, digits=2L, width,
                        xlab, slab, mlab, ilab, ilab.xpos, ilab.pos, order,
                        transf, atransf, targs, rows,
                        efac=1, pch=15, psize, col, border, lty, fonts,
                        cex, cex.lab, cex.axis, annosym, ...) {
  
  #########################################################################
  
  mstyle <- .get.mstyle("crayon" %in% .packages())
  
  if (!inherits(x, "rma"))
    stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))
  
  if (inherits(x, "rma.ls"))
    stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))
  
  na.act <- getOption("na.action")
  
  if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
    stop(mstyle$stop("Unknown 'na.action' specified under options()."))
  
  #if (!is.null(order))
  #   order <- match.arg(order, c("obs", "fit", "prec", "resid", "rstandard", "abs.resid", "abs.rstandard"))
  
  if (missing(transf))
    transf <- FALSE
  
  if (missing(atransf))
    atransf <- FALSE
  
  transf.char  <- deparse(substitute(transf))
  atransf.char <- deparse(substitute(atransf))
  
  if (is.function(transf) && is.function(atransf))
    stop(mstyle$stop("Use either 'transf' or 'atransf' to specify a transformation (not both)."))
  
  if (missing(targs))
    targs <- NULL
  
  if (missing(at))
    at <- NULL
  
  if (missing(ilab))
    ilab <- NULL
  
  if (missing(ilab.xpos))
    ilab.xpos <- NULL
  
  if (missing(ilab.pos))
    ilab.pos <- NULL
  
  if (missing(order))
    order <- NULL
  
  if (missing(psize))
    psize <- NULL
  
  if (missing(cex))
    cex <- NULL
  
  if (missing(cex.lab))
    cex.lab <- NULL
  
  if (missing(cex.axis))
    cex.axis <- NULL
  
  ### set default colors if user has not specified 'col' and 'border' arguments
  
  if (x$int.only) {
    
    if (missing(col)) {
      col <- c("black", "gray50") ### 1st color for summary polygon, 2nd color for credibility interval
    } else {
      if (length(col) == 1L)      ### if user only specified one value, assume it is for summary polygon
        col <- c(col, "gray50")
    }
    
    if (missing(border))
      border <- "black"           ### border color of summary polygon
    
  } else {
    
    if (missing(col))
      col <- "gray"               ### color of fitted values
    
    if (missing(border))
      border <- "gray"            ### border color of fitted values
    
  }
  
  if (missing(lty)) {
    lty <- c("solid", "dotted", "solid") ### 1st value = CIs, 2nd value = credibility interval, 3rd = horizontal line(s)
  } else {
    if (length(lty) == 1L)
      lty <- c(lty, "dotted", "solid")
    if (length(lty) == 2L)
      lty <- c(lty, "solid")
  }
  
  ### vertical expansion factor: 1st = CI end lines, 2nd = arrows, 3rd = summary polygon or fitted polygons
  
  if (length(efac) == 1L)
    efac <- rep(efac, 3)
  
  if (length(efac) == 2L)
    efac <- c(efac[1], efac[1], efac[2])
  
  ### annotation symbols vector
  
  if (missing(annosym))
    annosym <- c(" [", ", ", "]")
  if (length(annosym) != 3L)
    stop(mstyle$stop("Argument 'annosym' must be a vector of length 3."))
  
  level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))
  
  measure <- x$measure
  
  ### column header
  
  estlab <- .setlab(measure, transf.char, atransf.char, gentype=3, short=TRUE)
  if (is.expression(estlab)) {
    header.right <- parse(text=paste0("bold(", estlab, " * '", annosym[1], "' * '", 100*(1-level), "% CI'", " * '", annosym[3], "')"))
  } else {
    header.right <- paste0(estlab, annosym[1], 100*(1-level), "% CI", annosym[3])
  }
  
  if (is.logical(header)) {
    if (header) {
      header.left <- "Study"
    } else {
      header.left <- NULL
      header.right <- NULL
    }
  } else {
    if (is.expression(header)) {
      if (length(header) == 1L) {
        header.left <- parse(text=paste0("bold(", header, "')"))
      } else {
        header.left <- parse(text=paste0("bold(", header[1], ")"))
        header.right <- parse(text=paste0("bold(", header[2], ")"))
      }
    } else if (is.character(header)) {
      if (length(header) == 1L) {
        header.left <- header
      } else {
        header.left <- header[1]
        header.right <- header[2]
      }
    } else {
      stop(mstyle$stop("Argument 'header' must either be a logical or character vector."))
    }
    
  }
  
  ddd <- list(...)
  
  lplot     <- function(..., textpos) plot(...)
  labline   <- function(..., textpos) abline(...)
  lsegments <- function(..., textpos) segments(...)
  laxis     <- function(..., textpos) axis(...)
  lmtext    <- function(..., textpos) mtext(...)
  lpolygon  <- function(..., textpos) polygon(...)
  ltext     <- function(..., textpos) text(...)
  lpoints   <- function(..., textpos) points(...)
  
  ### TODO: remove this when there is a weights() function for 'rma.glmm' objects
  if (inherits(x, "rma.glmm") && showweights)
    stop(mstyle$stop("Option 'showweights=TRUE' not possible for 'rma.glmm' objects."))
  
  if (!is.null(ddd$subset))
    stop(mstyle$stop("Function does not have a 'subset' argument (could use 'order' argument instead)."))
  
  #########################################################################
  
  ### digits[1] for annotations, digits[2] for x-axis labels
  ### note: digits can also be a list (e.g., digits=list(2L,3))
  
  if (length(digits) == 1L)
    digits <- c(digits,digits)
  
  ### extract data and study labels
  ### note: yi.f/vi.f and pred may contain NAs
  
  yi <- x$yi.f
  vi <- x$vi.f
  X  <- x$X.f
  
  k <- length(yi)                              ### length of yi.f
  
  if (missing(slab)) {
    if (x$slab.null) {
      slab <- paste("Study", x$slab)         ### x$slab is always of length yi.f (i.e., NAs also have an slab)
    } else {
      slab <- x$slab                         ### note: slab must have same length as yi.f in rma object
    }                                         ### even when fewer studies used for model fitting (due to NAs)
  } else {
    if (length(slab) == 1L && is.na(slab))
      slab <- rep("", k)
  }
  
  if (length(yi) != length(slab))
    stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'slab' argument."))
  
  if (is.null(dim(ilab)))                      ### note: ilab must have same length as yi.f in rma object
    ilab <- cbind(ilab)                       ### even when fewer studies used for model fitting
  
  if (length(pch) == 1L)                       ### note: pch must have same length as yi.f in rma object
    pch <- rep(pch, k)                        ### or be equal to a single value (which is then repeated)
  
  if (length(pch) != length(yi))
    stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'pch' argument."))
  
  ### extract fitted values
  
  options(na.action = "na.pass")               ### using na.pass to get the entire vector (length of yi.f)
  
  if (x$int.only) {
    pred <- fitted(x)
    pred.ci.lb <- rep(NA_real_, k)
    pred.ci.ub <- rep(NA_real_, k)
  } else {
    temp <- predict(x, level=level)
    pred <- temp$pred
    if (addcred) {
      pred.ci.lb <- temp$cr.lb
      pred.ci.ub <- temp$cr.ub
    } else {
      pred.ci.lb <- temp$ci.lb
      pred.ci.ub <- temp$ci.ub
    }
  }
  
  if (inherits(x, "rma.glmm")) {            ### TODO: change this when there is a weights() function for 'rma.glmm' objects
    #weights <- NULL
    weights <- rep(1, k)
  } else {
    weights <- weights(x)                  ### these are the weights used for the actual model fitting
  }
  
  options(na.action = na.act)
  
  ### if user has set the point sizes
  
  if (!is.null(psize)) {                       ### note: psize must have same length as yi.f (including NAs)
    if (length(psize) == 1L)                  ### or be equal to a single value (which is then repeated)
      psize <- rep(psize, k)
    if (length(psize) != length(yi))
      stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'psize' argument."))
  }
  
  ### sort the data if requested
  
  if (!is.null(order)) {
    
    if (is.character(order)) {
      
      if (length(order) != 1L)
        stop(mstyle$stop("Incorrect length of 'order' argument."))
      
      if (order == "obs")
        sort.vec <- order(yi)
      if (order == "fit")
        sort.vec <- order(pred)
      if (order == "prec")
        sort.vec <- order(vi, yi)
      if (order == "resid")
        sort.vec <- order(yi-pred, yi)
      if (order == "rstandard")
        sort.vec <- order(rstandard(x)$z, yi)
      if (order == "abs.resid")
        sort.vec <- order(abs(yi-pred), yi)
      if (order == "abs.rstandard")
        sort.vec <- order(abs(rstandard(x)$z), yi)
      
    } else {
      sort.vec <- order                      ### in principle, can also subset with the order argument
    }
    
    yi         <- yi[sort.vec]
    vi         <- vi[sort.vec]
    X          <- X[sort.vec,,drop=FALSE]
    slab       <- slab[sort.vec]
    ilab       <- ilab[sort.vec,,drop=FALSE]  ### if ilab is still NULL, then this remains NULL
    pred       <- pred[sort.vec]
    pred.ci.lb <- pred.ci.lb[sort.vec]
    pred.ci.ub <- pred.ci.ub[sort.vec]
    weights    <- weights[sort.vec]
    pch        <- pch[sort.vec]
    psize      <- psize[sort.vec]             ### if psize is still NULL, then this remains NULL
    
  }
  
  k <- length(yi)                              ### in case length of k has changed
  
  ### set rows value
  
  if (missing(rows)) {
    rows <- k:1
  } else {
    if (length(rows) == 1L) {                 ### note: rows must be a single value or the same
      rows <- rows:(rows-k+1)                ### length of yi.f (including NAs) *after ordering/subsetting*
    }
  }
  
  if (length(rows) != length(yi))
    stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'rows' argument."))
  
  ### reverse order
  
  yi         <- yi[k:1]
  vi         <- vi[k:1]
  X          <- X[k:1,,drop=FALSE]
  slab       <- slab[k:1]
  ilab       <- ilab[k:1,,drop=FALSE]          ### if ilab is still NULL, then this remains NULL
  pred       <- pred[k:1]
  pred.ci.lb <- pred.ci.lb[k:1]
  pred.ci.ub <- pred.ci.ub[k:1]
  weights    <- weights[k:1]
  pch        <- pch[k:1]
  psize      <- psize[k:1]                     ### if psize is still NULL, then this remains NULL
  rows       <- rows[k:1]
  
  ### check for NAs in yi/vi and act accordingly
  
  yiviX.na <- is.na(yi) | is.na(vi) | apply(is.na(X), 1, any)
  
  if (any(yiviX.na)) {
    
    not.na <- !yiviX.na
    
    if (na.act == "na.omit") {
      yi         <- yi[not.na]
      vi         <- vi[not.na]
      X          <- X[not.na,,drop=FALSE]
      slab       <- slab[not.na]
      ilab       <- ilab[not.na,,drop=FALSE] ### if ilab is still NULL, then this remains NULL
      pred       <- pred[not.na]
      pred.ci.lb <- pred.ci.lb[not.na]
      pred.ci.ub <- pred.ci.ub[not.na]
      weights    <- weights[not.na]
      pch        <- pch[not.na]
      psize      <- psize[not.na]            ### if psize is still NULL, then this remains NULL
      
      rows.new <- rows                       ### rearrange rows due to NAs being omitted from plot
      rows.na  <- rows[!not.na]              ### shift higher rows down according to number of NAs omitted
      for (j in seq_len(length(rows.na))) {
        rows.new[rows >= rows.na[j]] <- rows.new[rows >= rows.na[j]] - 1
      }
      rows <- rows.new[not.na]
      
    }
    
    if (na.act == "na.fail")
      stop(mstyle$stop("Missing values in results."))
    
  }                                            ### note: yi/vi may be NA if na.act == "na.exclude" or "na.pass"
  
  k <- length(yi)                              ### in case length of k has changed
  
  ### calculate individual CI bounds
  
  ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
  ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
  
  ### if requested, apply transformation to yi's and CI bounds
  
  if (is.function(transf)) {
    if (is.null(targs)) {
      yi         <- sapply(yi, transf)
      ci.lb      <- sapply(ci.lb, transf)
      ci.ub      <- sapply(ci.ub, transf)
      pred       <- sapply(pred, transf)
      pred.ci.lb <- sapply(pred.ci.lb, transf)
      pred.ci.ub <- sapply(pred.ci.ub, transf)
    } else {
      yi         <- sapply(yi, transf, targs)
      ci.lb      <- sapply(ci.lb, transf, targs)
      ci.ub      <- sapply(ci.ub, transf, targs)
      pred       <- sapply(pred, transf, targs)
      pred.ci.lb <- sapply(pred.ci.lb, transf, targs)
      pred.ci.ub <- sapply(pred.ci.ub, transf, targs)
    }
  }
  
  ### make sure order of intervals is always increasing
  
  tmp <- .psort(ci.lb, ci.ub)
  ci.lb <- tmp[,1]
  ci.ub <- tmp[,2]
  
  tmp <- .psort(pred.ci.lb, pred.ci.ub)
  pred.ci.lb <- tmp[,1]
  pred.ci.ub <- tmp[,2]
  
  ### apply ci limits if specified
  
  if (!missing(clim)) {
    clim <- sort(clim)
    if (length(clim) != 2L)
      stop(mstyle$stop("Argument 'clim' must be of length 2."))
    ci.lb[ci.lb < clim[1]] <- clim[1]
    ci.ub[ci.ub > clim[2]] <- clim[2]
    pred.ci.lb[pred.ci.lb < clim[1]] <- clim[1]
    pred.ci.ub[pred.ci.ub > clim[2]] <- clim[2]
  }
  
  ### set default point sizes (if not specified by user)
  
  if (is.null(psize)) {
    # if (is.null(weights)) {
    #    if (any(vi <= 0, na.rm=TRUE)) {           ### in case any vi value is zero
    #       psize <- rep(1, k)
    #    } else {                                  ### default psize is proportional to inverse standard error
    #       wi    <- 1/sqrt(vi)                    ### note: vi's that are NA are ignored (but vi's whose yi is
    #       psize <- wi/sum(wi, na.rm=TRUE)        ### NA are NOT ignored; an unlikely case in practice)
    #       psize <- (psize - min(psize, na.rm=TRUE)) / (max(psize, na.rm=TRUE) - min(psize, na.rm=TRUE))
    #       psize <- (psize * 1.0) + 0.5           ### note: only vi's that are still in the subset are used for determining the default point sizes
    #       if (all(is.na(psize)))                 ### if k=1, then psize is NA, so catch this (and maybe some other problems)
    #          psize <- rep(1, k)
    #    }
    # } else {
    wi    <- weights
    psize <- wi/sum(wi, na.rm=TRUE)
    rng   <- max(psize, na.rm=TRUE) - min(psize, na.rm=TRUE)
    if (rng <= .Machine$double.eps^0.5) {
      psize <- rep(1, k)
    } else {
      psize <- (psize - min(psize, na.rm=TRUE)) / rng
      psize <- (psize * 1.0) + 0.5
    }
    if (all(is.na(psize)))
      psize <- rep(1, k)
    # }
  }
  
  #########################################################################
  
  ### total range of CI bounds
  
  rng <- max(ci.ub, na.rm=TRUE) - min(ci.lb, na.rm=TRUE)
  
  if (annotate) {
    if (showweights) {
      plot.multp.l <- 2.00
      plot.multp.r <- 2.00
    } else {
      plot.multp.l <- 1.20
      plot.multp.r <- 1.20
    }
  } else {
    plot.multp.l <- 1.20
    plot.multp.r <- 0.40
  }
  
  ### set plot limits
  
  if (missing(xlim)) {
    xlim <- c(min(ci.lb, na.rm=TRUE) - rng * plot.multp.l, max(ci.ub, na.rm=TRUE) + rng * plot.multp.r)
    xlim <- round(xlim, digits[[2]])
    #xlim[1] <- xlim[1]*max(1, digits[[2]]/2)
    #xlim[2] <- xlim[2]*max(1, digits[[2]]/2)
  }
  
  ### set x axis limits (at argument overrides alim argument)
  
  alim.spec <- TRUE
  
  if (missing(alim)) {
    if (is.null(at)) {
      alim <- range(pretty(x=c(min(ci.lb, na.rm=TRUE), max(ci.ub, na.rm=TRUE)), n=steps-1))
      alim.spec <- FALSE
    } else {
      alim <- range(at)
    }
  }
  
  ### make sure the plot and x axis limits are sorted
  
  alim <- sort(alim)
  xlim <- sort(xlim)
  
  ### plot limits must always encompass the yi values
  
  if (xlim[1] > min(yi, na.rm=TRUE)) { xlim[1] <- min(yi, na.rm=TRUE) }
  if (xlim[2] < max(yi, na.rm=TRUE)) { xlim[2] <- max(yi, na.rm=TRUE) }
  
  ### x axis limits must always encompass the yi values (no longer required)
  
  #if (alim[1] > min(yi, na.rm=TRUE)) { alim[1] <- min(yi, na.rm=TRUE) }
  #if (alim[2] < max(yi, na.rm=TRUE)) { alim[2] <- max(yi, na.rm=TRUE) }
  
  ### plot limits must always encompass the x axis limits
  
  if (alim[1] < xlim[1]) { xlim[1] <- alim[1] }
  if (alim[2] > xlim[2]) { xlim[2] <- alim[2] }
  
  ### allow adjustment of position of study labels and annotations via textpos argument
  
  if (is.null(ddd$textpos))
    ddd$textpos <- c(xlim[1], xlim[2])
  
  if (length(ddd$textpos) != 2L)
    stop(mstyle$stop("Argument 'textpos' must be of length 2."))
  
  if (is.na(ddd$textpos[1]))
    ddd$textpos[1] <- xlim[1]
  
  if (is.na(ddd$textpos[2]))
    ddd$textpos[2] <- xlim[2]
  
  ### set y axis limits
  
  if (missing(ylim)) {
    if (x$int.only && addfit) {
      ylim <- c(-1.5, k+top)
    } else {
      ylim <- c(0.5, k+top)
    }
  } else {
    ylim <- sort(ylim)
  }
  
  ### generate x axis positions if none are specified
  
  if (is.null(at)) {
    if (alim.spec) {
      at <- seq(from=alim[1], to=alim[2], length.out=steps)
    } else {
      at <- pretty(x=c(min(ci.lb, na.rm=TRUE), max(ci.ub, na.rm=TRUE)), n=steps-1)
    }
  } else {
    at[at < alim[1]] <- alim[1] ### remove at values that are below or above the axis limits
    at[at > alim[2]] <- alim[2]
    at <- unique(at)
  }
  
  ### x axis labels (apply transformation to axis labels if requested)
  
  at.lab <- at
  
  if (is.function(atransf)) {
    if (is.null(targs)) {
      at.lab <- formatC(sapply(at.lab, atransf), digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE))
    } else {
      at.lab <- formatC(sapply(at.lab, atransf, targs), digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE))
    }
  } else {
    at.lab <- formatC(at.lab, digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE))
  }
  
  #########################################################################
  
  ### set/get fonts (1st for study labels, 2nd for annotations, 3rd for ilab)
  ### when passing a named vector, the names are for 'family' and the values are for 'font'
  
  if (missing(fonts)) {
    fonts <- rep(par("family"), 3)
  } else {
    if (length(fonts) == 1L)
      fonts <- rep(fonts, 3)
    if (length(fonts) == 2L)
      fonts <- c(fonts, fonts[1])
  }
  
  if (is.null(names(fonts)))
    fonts <- structure(c(1L,1L,1L), names=fonts)
  
  par(family=names(fonts)[1], font=fonts[1])
  
  ### adjust margins
  
  par.mar <- par("mar")
  par.mar.adj <- par.mar - c(0,3,1,1)
  par.mar.adj[par.mar.adj < 0] <- 0
  par(mar = par.mar.adj)
  on.exit(par(mar = par.mar))
  
  ### start plot
  
  lplot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", bty="n", ...)
  
  ### horizontal title line
  
  labline(h=ylim[2]-(top-1), lty=lty[3], ...)
  
  ### add reference line
  
  if (is.numeric(refline))
    lsegments(refline, ylim[1]-5, refline, ylim[2]-(top-1), lty="dotted", ...)
  
  ### set cex, cex.lab, and cex.axis sizes as a function of the height of the figure
  
  par.usr <- par("usr")
  height  <- par.usr[4] - par.usr[3]
  
  if (is.null(cex)) {
    lheight <- strheight("O")
    cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * k * lheight), 1)
  }
  
  if (is.null(cex)) {
    cex <- par("cex") * cex.adj
  } else {
    if (is.null(cex.lab))
      cex.lab <- cex
    if (is.null(cex.axis))
      cex.axis <- cex
  }
  if (is.null(cex.lab))
    cex.lab <- par("cex") * cex.adj
  if (is.null(cex.axis))
    cex.axis <- par("cex") * cex.adj
  
  #########################################################################
  
  ### if addfit and not an intercept-only model, add fitted polygons
  
  if (addfit && !x$int.only) {
    
    for (i in seq_len(k)) {
      
      if (is.na(pred[i]))
        next
      
      lpolygon(x=c(max(pred.ci.lb[i], alim[1]), pred[i], min(pred.ci.ub[i], alim[2]), pred[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac[3], rows[i], rows[i]-(height/100)*cex*efac[3]), col=col, border=border, ...)
      
      ### this would only draw intervals if bounds fall within alim range
      #if ((pred.ci.lb[i] > alim[1]) && (pred.ci.ub[i] < alim[2]))
      #   lpolygon(x=c(pred.ci.lb[i], pred[i], pred.ci.ub[i], pred[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac[3], rows[i], rows[i]-(height/100)*cex*efac[3]), col=col, border=border, ...)
      
    }
    
  }
  
  #########################################################################
  
  ### if addfit and intercept-only model, add fixed/random-effects model polygon
  
  if (addfit && x$int.only) {
    
    if (inherits(x, "rma.mv") && x$withG && x$tau2s > 1) {
      
      if (!is.logical(addcred)) {
        ### for multiple tau^2 (and gamma^2) values, need to specify level(s) of the inner factor(s) to compute the credibility interval
        ### this can be done via the addcred argument (i.e., instead of using a logical, one specifies the level(s))
        if (length(addcred) == 1L)
          addcred <- c(addcred, addcred)
        temp <- predict(x, level=level, tau2.levels=addcred[1], gamma2.levels=addcred[2])
        addcred <- TRUE ### set addcred to TRUE, so if (x$method != "FE" && addcred) further below works
      } else {
        if (addcred) {
          ### here addcred=TRUE, but user has not specified the level, so throw an error
          stop(mstyle$stop("Need to specify the level of the inner factor(s) via the 'addcred' argument."))
        } else {
          ### here addcred=FALSE, so just use the first tau^2 and gamma^2 arbitrarily (so predict() works)
          temp <- predict(x, level=level, tau2.levels=1, gamma2.levels=1)
        }
      }
      
    } else {
      
      temp <- predict(x, level=level)
      
    }
    
    beta       <- temp$pred
    beta.ci.lb <- temp$ci.lb
    beta.ci.ub <- temp$ci.ub
    beta.cr.lb <- temp$cr.lb
    beta.cr.ub <- temp$cr.ub
    
    if (is.function(transf)) {
      if (is.null(targs)) {
        beta       <- sapply(beta, transf)
        beta.ci.lb <- sapply(beta.ci.lb, transf)
        beta.ci.ub <- sapply(beta.ci.ub, transf)
        beta.cr.lb <- sapply(beta.cr.lb, transf)
        beta.cr.ub <- sapply(beta.cr.ub, transf)
      } else {
        beta       <- sapply(beta, transf, targs)
        beta.ci.lb <- sapply(beta.ci.lb, transf, targs)
        beta.ci.ub <- sapply(beta.ci.ub, transf, targs)
        beta.cr.lb <- sapply(beta.cr.lb, transf, targs)
        beta.cr.ub <- sapply(beta.cr.ub, transf, targs)
      }
    }
    
    ### make sure order of intervals is always increasing
    
    tmp <- .psort(beta.ci.lb, beta.ci.ub)
    beta.ci.lb <- tmp[,1]
    beta.ci.ub <- tmp[,2]
    
    tmp <- .psort(beta.cr.lb, beta.cr.ub)
    beta.cr.lb <- tmp[,1]
    beta.cr.ub <- tmp[,2]
    
    ### apply ci limits if specified
    
    if (!missing(clim)) {
      beta.ci.lb[beta.ci.lb < clim[1]] <- clim[1]
      beta.ci.ub[beta.ci.ub > clim[2]] <- clim[2]
      beta.cr.lb[beta.cr.lb < clim[1]] <- clim[1]
      beta.cr.ub[beta.cr.ub > clim[2]] <- clim[2]
    }
    
    ### add credibility interval
    
    if (x$method != "FE" && addcred) {
      
      lsegments(max(beta.cr.lb, alim[1]), -1, min(beta.cr.ub, alim[2]), -1, lty=lty[2], col=col[2], ...)
      
      if (beta.cr.lb >= alim[1]) {
        lsegments(beta.cr.lb, -1-(height/150)*cex*efac[1], beta.cr.lb, -1+(height/150)*cex*efac[1], col=col[2], ...)
      } else {
        lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(-1, -1+(height/150)*cex*efac[2], -1-(height/150)*cex*efac[2], -1), col=col[2], border=col[2], ...)
      }
      
      if (beta.cr.ub <= alim[2]) {
        lsegments(beta.cr.ub, -1-(height/150)*cex*efac[1], beta.cr.ub, -1+(height/150)*cex*efac[1], col=col[2], ...)
      } else {
        lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(-1, -1+(height/150)*cex*efac[2], -1-(height/150)*cex*efac[2], -1), col=col[2], border=col[2], ...)
      }
      
    }
    
    ### polygon for the summary estimate
    
    lpolygon(x=c(beta.ci.lb, beta, beta.ci.ub, beta), y=c(-1, -1+(height/100)*cex*efac[3], -1, -1-(height/100)*cex*efac[3]), col=col[1], border=border, ...)
    
    ### add label for model estimate
    
    if (missing(mlab))
      mlab <- ifelse((x$method=="FE"), "FE Model", "RE Model")
    
    ltext(ddd$textpos[1], -1, mlab, pos=4, cex=cex, ...)
    
  }
  
  #########################################################################
  
  ### add x axis
  
  laxis(side=1, at=at, labels=at.lab, cex.axis=cex.axis, ...)
  
  ### add x axis label
  
  if (missing(xlab))
    xlab <- .setlab(measure, transf.char, atransf.char, gentype=1)
  
  lmtext(xlab, side=1, at=min(at) + (max(at)-min(at))/2, line=par("mgp")[1]-0.5, cex=cex.lab, ...)
  
  ### add CI ends (either | or <> if outside of axis limits)
  
  for (i in seq_len(k)) {
    
    ### need to skip missings, as if() check below will otherwise throw an error
    if (is.na(yi[i]) || is.na(vi[i]))
      next
    
    ### if the lower bound is actually larger than upper x-axis limit, then everything is to the right and just draw a polygon pointing in that direction
    if (ci.lb[i] >= alim[2]) {
      lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col="black", ...)
      next
    }
    
    ### if the upper bound is actually lower than lower x-axis limit, then everything is to the left and just draw a polygon pointing in that direction
    if (ci.ub[i] <= alim[1]) {
      lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col="black", ...)
      next
    }
    
    lsegments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], alim[2]), rows[i], lty=lty[1], ...)
    
    if (ci.lb[i] >= alim[1]) {
      lsegments(ci.lb[i], rows[i]-(height/150)*cex*efac[1], ci.lb[i], rows[i]+(height/150)*cex*efac[1], ...)
    } else {
      lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col="black", ...)
    }
    
    if (ci.ub[i] <= alim[2]) {
      lsegments(ci.ub[i], rows[i]-(height/150)*cex*efac[1], ci.ub[i], rows[i]+(height/150)*cex*efac[1], ...)
    } else {
      lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col="black", ...)
    }
    
  }
  
  ### add study labels on the left
  
  ltext(ddd$textpos[1], rows, slab, pos=4, cex=cex, ...)
  
  ### add info labels
  
  if (!is.null(ilab)) {
    if (is.null(ilab.xpos))
      stop(mstyle$stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'."))
    if (length(ilab.xpos) != ncol(ilab))
      stop(mstyle$stop(paste0("Number of 'ilab' columns (", ncol(ilab), ") does not match length of 'ilab.xpos' argument (", length(ilab.xpos), ").")))
    if (!is.null(ilab.pos) && length(ilab.pos) == 1L)
      ilab.pos <- rep(ilab.pos, ncol(ilab))
    par(family=names(fonts)[3], font=fonts[3])
    for (l in seq_len(ncol(ilab))) {
      ltext(ilab.xpos[l], rows, ilab[,l], pos=ilab.pos[l], cex=cex, ...)
    }
    par(family=names(fonts)[1], font=fonts[1])
  }
  
  ### add study annotations on the right: yi [LB, UB]
  ### and add model fit annotations if requested: b [LB, UB]
  ### (have to add this here, so that alignment is correct)
  
  if (annotate) {
    
    if (is.function(atransf)) {
      
      if (is.null(targs)) {
        if (addfit && x$int.only) {
          annotext <- cbind(sapply(c(yi, beta), atransf), sapply(c(ci.lb, beta.ci.lb), atransf), sapply(c(ci.ub, beta.ci.ub), atransf))
        } else {
          annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, atransf), sapply(ci.ub, atransf))
        }
      } else {
        if (addfit && x$int.only) {
          annotext <- cbind(sapply(c(yi, beta), atransf, targs), sapply(c(ci.lb, beta.ci.lb), atransf, targs), sapply(c(ci.ub, beta.ci.ub), atransf, targs))
        } else {
          annotext <- cbind(sapply(yi, atransf, targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, atransf, targs))
        }
      }
      
      ### make sure order of intervals is always increasing
      
      tmp <- .psort(annotext[,2:3])
      annotext[,2:3] <- tmp
      
    } else {
      
      if (addfit && x$int.only) {
        annotext <- cbind(c(yi, beta), c(ci.lb, beta.ci.lb), c(ci.ub, beta.ci.ub))
      } else {
        annotext <- cbind(yi, ci.lb, ci.ub)
      }
      
    }
    
    if (showweights) {
      if (addfit && x$int.only) {
        annotext <- cbind(c(unname(weights),100), annotext)
      } else {
        annotext <- cbind(unname(weights), annotext)
      }
    }
    
    annotext <- formatC(annotext, format="f", digits=digits[[1]])
    
    if (missing(width)) {
      width <- apply(annotext, 2, function(x) max(nchar(x)))
    } else {
      if (length(width) == 1L)
        width <- rep(width, ncol(annotext))
    }
    
    for (j in seq_len(ncol(annotext))) {
      annotext[,j] <- formatC(annotext[,j], width=width[j])
    }
    
    if (showweights) {
      annotext <- cbind(annotext[,1], "%   ", annotext[,2], annosym[1], annotext[,3], annosym[2], annotext[,4], annosym[3])
    } else {
      annotext <- cbind(annotext[,1], annosym[1], annotext[,2], annosym[2], annotext[,3], annosym[3])
    }
    
    annotext <- apply(annotext, 1, paste, collapse="")
    annotext[grepl("NA", annotext, fixed=TRUE)] <- ""
    
    par(family=names(fonts)[2], font=fonts[2])
    if (addfit && x$int.only) {
      ltext(ddd$textpos[2], c(rows,-1), labels=annotext, pos=2, cex=cex, ...)
    } else {
      ltext(ddd$textpos[2], rows, labels=annotext, pos=2, cex=cex, ...)
    }
    par(family=names(fonts)[1], font=fonts[1])
    
  }
  
  ### add yi points
  
  for (i in seq_len(k)) {
    
    ### need to skip missings, as if() check below will otherwise throw an error
    if (is.na(yi[i]))
      next
    
    if (yi[i] >= alim[1] && yi[i] <= alim[2])
      lpoints(yi[i], rows[i], pch=pch[i], cex=cex*psize[i], ...)
    
  }
  
  #lpoints(yi, rows, pch=pch, cex=cex*psize, ...)
  
  ### add horizontal line at 0 for the standard FE/RE model display
  
  if (x$int.only && addfit)
    labline(h=0, lty=lty[3], ...)
  
  ### add header
  
  ltext(ddd$textpos[1], ylim[2]-(top-1)+1, header.left, pos=4, font=2, cex=cex, ...)
  ltext(ddd$textpos[2], ylim[2]-(top-1)+1, header.right, pos=2, font=2, cex=cex, ...)
  
  #########################################################################
  
  ### return some information about plot invisibly
  
  res <- list('xlim'=par("usr")[1:2], 'alim'=alim, 'at'=at, 'ylim'=ylim, 'rows'=rows, 'cex'=cex, 'cex.lab'=cex.lab, 'cex.axis'=cex.axis)
  
  invisible(res)
  
}