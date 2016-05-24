#------------------------------------------------------------------------------
# function to calculate the CV from a given 90% confidence interval
# 
# Author: dlabes, hschuetz
# adapted to unbalanced studies by Benjamin Lang
#------------------------------------------------------------------------------

CVfromCI <- function(point, lower, upper, n, design = "2x2", alpha = 0.05, 
                     robust = FALSE) 
{
  if (missing(n)) stop("Sample size n must be given!", call. = FALSE)
  # according to Helmut's suggestion
  if ((missing(lower) && missing(upper)) || ((missing(lower) && 
       missing(point)) || (missing(upper) && missing(point)))) {
    stop("At least both CLs or PE and one CL must be given!", call. = FALSE)
  }
  # calculate missing lower or upper
  if (!missing(point) && (missing(lower) || missing(upper))) {
    ifelse(missing(lower), lower <- point^2/upper, upper <- point^2/lower)
  }
  # calculate missing PE
  if (missing(point)) point <- sqrt(lower * upper)
  # should we really do that?
  if (length(point) > 1 || length(lower) > 1 || length(upper) > 1)
    stop("point, lower, upper must be scalars.", call. = FALSE)
  # check design
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ", design, " unknown!", call. = FALSE)
  ades <- .design.props(d.no)
  if (length(n) == 1) {
    # n given as ntotal
    n <- nvec(n = n, grps = ades$steps)
    if (n[1] != n[length(n)]) {
      message("Unbalanced ", design, " design. n(i)= ", paste(n, collapse = "/"),
              " assumed.")
    }
  }
  else {
    # n given as vector of # of subjects in (sequence) groups    
    if (length(n) != ades$steps) stop("Length of n vector must be ", ades$steps, "!")
  }
  nc <- sum(1/n)
  n <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  df <- eval(.design.df(ades, robust = robust))
  tval <- qt(1 - alpha, df)
  s1 <- (log(point) - log(lower))/se.fac/tval
  s2 <- (log(upper) - log(point))/se.fac/tval
  sw <- 0.5 * (s1 + s2)
  if (abs(s1 - s2)/sw > 0.1) {
    warning(paste("sigma based on pe & lower CL more than 10% different than\n", 
            "sigma based on pe & upper CL. Check input."), call. = FALSE)
  }
  return(se2CV(sw))
}
# ---------------------------------------------------------------------------
# alias to CVfromCI
CI2CV <- function(point, lower, upper, n, design="2x2", alpha=0.05, 
                  robust=FALSE)
{
  CVfromCI(point, lower, upper, n, design=design, alpha=alpha, 
           robust=robust)
}
