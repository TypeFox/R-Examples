# ----------------------------------------------------------------------------
# calculate the TOST p-value(s)
# contributed by Benjamin Lang
# reviewed by D.Labes
# ----------------------------------------------------------------------------
pvalues.TOST <- function(pe, CV, n, logscale = TRUE, theta1, theta2,
                        design = "2x2", robust = FALSE, both=TRUE)
{
  pvalue.TOST(pe=pe,CV=CV,n=n,logscale=logscale, theta1=theta1, theta2=theta2,
              design = design, robust = robust, both=both)
}
  
pvalue.TOST <- function(pe, CV, n, logscale = TRUE, theta1, theta2,
                        design = "2x2", robust = FALSE, both=FALSE)
{
  # Calculates the p-value(s) for the TOST procedure
  #
  # Args:
  #   pe:     Observed point estimate (ratio of geometric means)
  #   CV:     Observed coefficient of variation as ratio
  #   n:      Number of subjects
  #   theta1: Lower bioequivalence limit
  #   theta2: Upper bioequivalence limit
  #   design: Character string describing the study design
  #   robust: Indicates whether robust degrees of freedem should be used
  #   both:   Indicates if both p-values or the maximum of them will be given back
  #
  # Returns:
  #   p-value(s)
  
  # Error handling
  # these checks are not really necessary, try to use without
  if (missing(pe)) stop("Point estimator pe must be given!")
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")
  
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ", design, " unknown!", call. = FALSE)
  ades <- .design.props(d.no)
  dfe  <- .design.df(ades, robust = robust)
  # we use always bkni
  #bk   <- ades$bk
  if (length(n) == 1) {
    # total n given, but may lead to unbalanced design    
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance (function nvec() from Helper_dp.R)
    n <- nvec(n=n, grps=ades$steps)
    if (n[1]!=n[length(n)]){
      message("Unbalanced ",design, " design. n(i)= ", paste(n, collapse="/"), 
              " assumed.")
    } 
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
  }
  nc <- sum(1/n)
  n  <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  
  df <- eval(dfe)
  if (any(df < 1)) stop("n(s) too small, degrees of freedom < 1!")
  
  if (logscale) {
    if (any(pe <= 0))  stop("pe(s) must be > 0!")
    if (missing(theta1)) theta1 <- 0.8
    if (missing(theta2)) theta2 <- 1/theta1
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff   <- log(pe)
    if (any(CV <= 0)) stop("CV must be >= 0!")
    # se of ldiff (sem)
    se <- CV2se(CV) * se.fac
  } else {
    # cat() changed to message to allow suppressMessages()
    message("It is assumed that CV is the standard deviation.")
    if (missing(theta1)) theta1 <- log(0.8)
    if (missing(theta2)) theta2 <- -theta1
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff   <- pe
    se <- CV * se.fac
  }
  # Calculate value of test statistics
  t1 <- (ldiff - ltheta1) / se
  t2 <- (ldiff - ltheta2) / se
  # Calculate P(T > t1 | H0) and P(T < t2 | H0)
  p.left  <- pt(t1, df = df, lower.tail = FALSE)  # right-tailed
  p.right <- pt(t2, df = df, lower.tail = TRUE)   # left-tailed
  ret <- cbind(p.left, p.right)    # gives a matrix if p.left/p.right are vectors
  #colnames(ret) <- c("p(\u2264\u03981)","p(\u2265\u03982)")
  if (nrow(ret)==1) ret <- ret[1,] # gives a vector
  # Both null hypotheses need to be rejected simultaneously
  # thus Benjamin prefers max of both pvalues
  if (both) return(ret) else return(pmax(p.left, p.right))
}
