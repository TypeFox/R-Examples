ns.test <-
function(x, g = NULL, q, alternative = c("two.sided","less","greater"),
                    paired = FALSE, compared = FALSE, alpha = 0.05)
{
  alternative <- match.arg(alternative)
  if (!is.numeric(x))
    stop("'x' must be numeric") 
  if (is.null(g))
  {
    DNAME <- deparse(substitute(x))
    METHOD <- "One-sample normal score test"
    if (missing(q)) stop("please specify 'q'")
    if (!missing(q) && ((length(q) > 1L) || !is.finite(q))) 
      stop("'q' must be a finite single number")
    if (paired) 
      stop("'g' is missing for paired test")
    x <- x[is.finite(x)]
  }
  else
  {
    gr <- length(unique(g))
    if (length(g) < 2 || gr < 2)
      stop(paste("'g' must have at least two groups or please specify q = ", g, 
                 "for one-sample test"))                
    else if (gr == 2)
    {
      DNAME <- paste(deparse(substitute(x)), "with group", deparse(substitute(g)))
      METHOD <- "Two-sample normal score test"
      if (!paired && !missing(q))
        stop("argument 'q' can only be used for one-sample test")
      xx <- x[g == unique(g)[1]]; yy <- x[g == unique(g)[2]]
      x <- NULL; y <- NULL; x <- xx[is.finite(xx)]; y <- yy[is.finite(yy)]
      n <- as.integer(length(y))
      if (paired)
      {
        METHOD <- "Paired two-sample normal score test"
        if (length(x) != n)
          stop("two groups must have the same length of observations")
        OK <- complete.cases(x,y)
        x <- x[OK] - y[OK]
        y <- NULL  
      }
    }
    else 
    {
      DNAME <- paste(deparse(substitute(x)), "with group", deparse(substitute(g)))
      METHOD <- "Multiple-sample normal score test"
      if (!missing(q))
        stop("argument 'q' can only be used for one-sample test")
      if (paired)
        stop("argument 'paired' can only be used for two-sample test")
    }    
  }
  if (is.null(g) || (!is.null(g) && length(unique(g)) < 3))
  {
    if (compared)
      stop("argument 'compared' can be only used for multiple-sample test")
    if (!missing(alpha))
      stop("argument 'alpha' can be only used for multiple-sample test")
  }
  m <- as.integer(length(x))
  if (m < 1L) 
    stop("not enough 'x' data")
  if (is.null(g) || (!is.null(g) && paired))
  {
    if (!is.null(g) && paired && !missing(q)) 
      warning("'q' has been setted to be 0 for paired two-sample test")
    q <- ifelse((!is.null(g) && paired), 0, q)
    x <- x - q
    zeros <- any(x == 0)
    if (zeros) 
      x <- x[x != 0]
    m <- as.integer(length(x))
    if (m < 1L) 
      stop("not enough data for first group")
    s <- qnorm((1 + rank(abs(x))/(m + 1))/2)*sign(x)
    STAT <- sum(s)/sqrt(sum(s^2)) 
    P.VAL <- switch(alternative, two.sided = 2*min(pnorm(STAT), 1 - pnorm(STAT)),
                    less = pnorm(STAT), greater = 1 - pnorm(STAT))
    ALTERNATIVE <- switch(alternative, two.sided = paste("true location is not equal to", 
                          q),less = paste("true location is less than", q), 
                          greater = paste("true location is greater than", q))
    RESULT <- list(p.value = P.VAL, alternative = ALTERNATIVE)
  }
  else if (gr == 2)
  {
    if (n < 1L) 
      stop("not enough data for second group")
    w <- qnorm(rank(c(x,y))/(m + n + 1))
    r.x <- w[1:m]; r.y <- w[(m + 1):(m + n)]
    STAT <- sum(w[1:m])/sqrt(m*n*sum(w^2)/((m + n - 1)*(m + n)))
    P.VAL <- switch(alternative, two.sided = 2*min(pnorm(STAT), 1 - pnorm(STAT)),
                    less = pnorm(STAT), greater = 1 - pnorm(STAT))
    ALTERNATIVE <- switch(alternative,two.sided = "true location shift is not equal to 0",
                          less = "true location shift is less than 0",
                          greater = "true location shift is greater than 0")
    RESULT <- list(p.value = P.VAL, alternative = ALTERNATIVE)
  }
  else
  {
    group <- unique(g)
    w <- qnorm(rank(x)/(m + 1))
    numerator <- 0; denominator <- 0
    for (i in 1:gr)
    {
      numerator <- numerator + (sum(w[g == group[i]]))^2/length(w[g == group[i]])
      denominator <- denominator + sum((w[g == group[i]])^2)
    }
    STAT <- (m - 1)*numerator/denominator
    P.VAL <- 1 - pchisq(STAT, gr - 1)
    ALTERNATIVE <- "at least one population is greater than at least one of the others"
    PARAMETER <- setNames(gr - 1, "df")
    RESULT <- list(p.value = P.VAL, alternative = ALTERNATIVE, parameter = PARAMETER)
    if (compared)
    {
      if (alpha < 0 || alpha > 1)
        stop("significant level alpha must be between 0 and 1")
      A <- NULL; n <- NULL
      for (i in 1:gr) 
      {
        n[i] <- length(w[g == group[i]])
        A[i] <- mean(w[g == group[i]])
      }  
      COMP <- matrix(NA,gr,gr)
      for (i in 1:(gr - 1))
        for (j in (i + 1):gr)
        {
          COMP[i,j] <- ifelse(abs(A[i] - A[j]) > sqrt(denominator/(m - 1))*
                                sqrt((m - 1 - STAT)/(m - gr))*sqrt(1/n[i] + 1/n[j])*
                                qt(1 - alpha/2/gr, n[i] + n[j] - 2), "YES","NO")
        }
      INDEX.YES <- which(COMP == "YES", arr.ind = TRUE)
      INDEX.NO <- which(COMP == "NO", arr.ind = TRUE)
      COMPARE <- rbind(cbind(INDEX.YES,matrix(rep("YES",nrow(INDEX.YES)))),
                       cbind(INDEX.NO,matrix(rep("NO",nrow(INDEX.NO)))))
      colnames(COMPARE) <- c("group", "group", "DIFF")
      RESULT <- c(RESULT, list(compare = COMPARE))
    }
  }
  STATISTIC <- setNames(STAT, "statistic")  
  RVAL <- c(list(data.name = DNAME, method = METHOD,
                 statistic = STATISTIC), RESULT)
  if (!compared) class(RVAL) <- "htest"
  return(RVAL)
}
