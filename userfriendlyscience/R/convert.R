### See http://statistiki.eu/wiki/Converting_Effect_Sizes for some formulae

###########################################################################
### Converting from: Pearson's r
###########################################################################

convert.r.to.t <- function(r, n) {
  return(r * sqrt((n - 2) / (1-r^2)));
}

convert.r.to.d <- function(r) {
  return( (r*2) / sqrt(1 - r^2));
}

convert.r.to.p <- function(r, n) {
  t <- convert.r.to.t(r, n);
  return(convert.t.to.p(t, n - 2));
}

convert.r.to.fisherz <- function(r) {
  return(.5 * log((1+r) / (1-r)));
}

###########################################################################
### Converting from: Student's t
###########################################################################

convert.t.to.r <- function(t, n) {
  return(t / (sqrt(n-2+t^2)));
}

convert.t.to.p <- function(t, df) {
  return(2*pt(-abs(t),df));
}

convert.t.to.d <- function(t, df=NULL, n1=NULL, n2=NULL, proportion=.5) {
 
  if (is.null(df) && !is.null(n1) && !is.null(n2)) {
    groupSize1 <- n1;
    groupSize2 <- n2;
  }
  else if (!is.null(df) && is.null(n1) && is.null(n2)) {
    groupSize1 <-      proportion  * (df + 2);
    groupSize2 <- (1 - proportion) * (df + 2);
  }
  else {
    warning("Specify either df (and ideally proportion) or n1 and n2! Returning NA.");
    return(NA);
  }
  
  d <- t * sqrt(((groupSize1 + groupSize2) / (groupSize1 * groupSize2)) *
                ((groupSize1 + groupSize2) / (groupSize1 + groupSize2 - 2)));
  
  return(d);
}

###########################################################################
### Converting from: Cohen's d
###########################################################################

convert.d.to.r <- function(d, n1 = NULL, n2 = NULL) {
  if (is.null(n1) && is.null(n2)) {
    a <- 4;
  } else {
    a <- (n1 + n2) ^ 2 / (n1 * n2)
  }
  return(d / sqrt(d^2 + a));
}

convert.d.to.t <- function(d, n) {
  return(sqrt(sqrt(n) * d));
}

convert.d.to.logodds <- function(d) {
  if (!is.numeric(d) || (length(d) > 1)) {
    stop("The 'd' argument is not a single numeric value!");
  }
  return(d * (pi / sqrt(3)));
}

convert.d.to.variance <- function(d, n1, n2) {
  return( (((n1+n2) / (n1*n2)) + ((d^2) / (2*(n1+n2-2)))) * ((n1+n2) / (n1+n2-2)) );
}

###########################################################################
### Converting from: Log Odds
###########################################################################

convert.logodds.to.d <- function(logodds) {
  return(logodds * (sqrt(3) / pi));
}

convert.logodds.to.r <- function(logodds) {
  return(convert.d.to.r(convert.logodds.to.d(logodds)));
}

###########################################################################
### Converting from: Odds Ratio
###########################################################################

convert.or.to.d <- function(or) {
  return(log(or) * (sqrt(3) / pi));
}

convert.or.to.r <- function(or) {
  return(convert.d.to.r(convert.logodds.to.d(log(or))));
}

###########################################################################
### Converting from: Chi Square
###########################################################################

convert.chisq.to.V <- function(chisq, n, minDim) {
  if (!is.numeric(chisq) || (length(chisq) > 1)) {
    stop("The 'chisq' argument is not a single numeric value!");
  }
  if (!is.numeric(n) || (length(n) > 1)) {
    stop("The 'n' argument is not a single numeric value!");
  }
  if (!is.numeric(minDim) || (length(minDim) > 1)) {
    stop("The 'minDim' argument is not a single numeric value!");
  }
  res <- as.numeric(sqrt(chisq/(n*(minDim - 1))));
  return(ifelse(is.finite(res), res, NA));
}

convert.chisq.to.p <- function(chisq, df, lower.tail=FALSE) {
  return(2*pchisq(chisq, df, lower.tail=lower.tail));
}

###########################################################################
### Converting from: F
###########################################################################

convert.f.to.p <- function(f, df1, df2, lower.tail=FALSE) {
  return(2*pf(f, df1, df2, lower.tail=lower.tail));
}

convert.f.to.d <- function(f, df1, df2 = NULL, n1=NULL, n2=NULL, proportion=.5) {
  if (df1 != 1) {
    warning("You can only convert an F value for the comparison of two groups to Cohen's d, ",
            "and you specified a df1 of ", df1, ", which means this F value concerns the comparison ",
            "of ", df1 + 1, " groups. Returning NA.");
    return(NA);
  }
  else if (is.null(df2) && !is.null(n1) && !is.null(n2)) {
    groupSize1 <- n1;
    groupSize2 <- n2;
  }
  else if (!is.null(df2) && is.null(n1) && is.null(n2)) {
    groupSize1 <- proportion * (df1 + df2 + 1);
    groupSize2 <- (1 - proportion) * (df1 + df2 + 1);
  }
  else {
    warning("Specify either df2 (and ideally proportion) or n1 and n2! Returning NA.");
    return(NA);
  }
  
  d <- sqrt(f * ((groupSize1 + groupSize2) / (groupSize1 * groupSize2)) *
              ((groupSize1 + groupSize2) / (groupSize1 + groupSize2 - 2)));
  
  return(d);
}

convert.f.to.etasq <- function(f, df1, df2) {
  return( (f * df1) / ((f * df1) + df2) );
}

convert.f.to.omegasq <- function(f, df1, df2) {
  return( (f - 1) / (f + (df2 + 1) / (df1)) );
}

###########################################################################
### Converting from: Eta square
###########################################################################

convert.etasq.to.cohensf <- function(etasq) {
  return(sqrt(etasq / (1-etasq)));
}

###########################################################################
### Converting from: Beta (regression weight)
###########################################################################

convert.b.to.t <- function(b, se) {
  return(b/se);
}

###########################################################################
### Converting from: Fisher's z
###########################################################################

convert.fisherz.to.r <- function(z) {
  return((exp(2 * z) - 1) / (exp(2*z)+1));
}
