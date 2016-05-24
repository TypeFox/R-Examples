"or.midp" <-
function(x, conf.level = 0.95, byrow = TRUE, interval = c(0, 1000)){
  ##housekeeping
  if(is.vector(x)){
    if(!is.numeric(x)){stop("vector must be numeric")}
    if(length(x)!=4){stop("vector must be of length 4")}
    x <- matrix(x, 2, 2, byrow = byrow)
  }
  if(is.matrix(x)){
    if(!is.numeric(x)){stop("matrix must be numeric")}
    if(nrow(x)!=2 || ncol(x)!=2){stop("must be 2 x 2 matrix")}
    a1 <- x[1,1]; a0 <- x[1,2]; b1 <- x[2,1]; b0 <- x[2,2]
  } else {stop("must be numeric vector of length=4 or 2x2 numeric matrix")}
  ##median-unbiased estimate function
  mue <- function(a1, a0, b1, b0, or){
    mm <- matrix(c(a1,a0,b1,b0),2,2, byrow=TRUE)
    fisher.test(mm, or=or, alternative="l")$p-fisher.test(x=x, or=or,
                             alternative="g")$p
  }
  ##mid-p function
  midp <- function(a1, a0, b1, b0, or = 1){
    mm <- matrix(c(a1,a0,b1,b0),2,2, byrow=TRUE)
    lteqtoa1 <- fisher.test(mm,or=or,alternative="l")$p.val
    gteqtoa1 <- fisher.test(mm,or=or,alternative="g")$p.val
    0.5*(lteqtoa1-gteqtoa1+1)
  }
  alpha <- 1 - conf.level
  ##root finding
  EST <- uniroot(function(or){
    mue(a1, a0, b1, b0, or)
  }, interval = interval)$root
  LCL <- uniroot(function(or){
    1-midp(a1, a0, b1, b0, or)-alpha/2
  },  interval = interval)$root
  UCL <- 1/uniroot(function(or){
    midp(a1, a0, b1, b0, or=1/or)-alpha/2
  },  interval = interval)$root
  rr <- list(x = x,
             estimate = EST,
             conf.int = c(LCL, UCL),
             conf.level = conf.level)
  attr(rr, "method") <- "median-unbiased estimate & mid-p exact CI"
  return(rr)
}
