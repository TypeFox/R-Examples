fp.test <- function(x,...) {
  UseMethod("fp.test")
}

fp.test.default <- function (x,y,delta=0,alternative="two.sided",...) {
  dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
  names(delta) <- "difference in location"
  y <- y - delta
  ds <- matrix(outer(y,x,"-"),ncol=1)
  m <- length(x)
  n <- length(y)
  mu <- m*n/2
  ind <- rep(0,m*n)
  ind[ds>0] <- 1
  ind[ds==0] <- 0.5
  ts <- sum(ind)-mu
  place <- function (x,y) {
    ic <- 0
    n <- length(y)
    ac <- 0
    ys <- sort(y)
    i <- 1
    while (ac==0) {
	tmp <- ys[i]
	if (x>tmp) {
	  i <- i+1
	  ic <- ic+1
	} else {
	  ac <- 1
	}
	if (i>n) {ac <- 1}
    }
    ic
  }
  pl <- apply(as.matrix(x),1,place,y)
  ql <- apply(as.matrix(y),1,place,x)
  v1 <- (m-1)*var(pl)
  v2 <- (n-1)*var(ql)
  sig <- sqrt(v1+v2+(mean(pl)*mean(ql)))
  zp <- ts/sig
  std <- zp
  names(std) <- "U*"
  zp <- (ts-0.5)/sig
  zn <- (ts+0.5)/sig
  if (alternative=="greater") {
    pval <- 1-pnorm(zp)
    zs <- zp
  } else if (alternative=="less") {
    pval <- pnorm(zn)
    zs <- zn
  } else if (alternative=="two.sided") {
    if (ts>=0) {
	pval <- 2*(1-pnorm(abs(zp)))
	zs <- zp
    } else {
	pval <- 2*pnorm(zn)
	zs <- zn
    }
  }
  result <- list(statistic=std,p.value=pval,alternative=alternative,method="Fligner-Policello test",
    data.name=dname,null.value=delta)
  class(result) <- "htest"
  return(result)
}

fp.test.formula <- function(formula,data,subset,...) {
  if (missing(formula) || (length(formula)!=3L) || (length(attr(terms(formula[-2L]), 
      "term.labels"))!=1L)) {stop("'formula' missing or incorrect")}
  m <- match.call(expand.dots=FALSE)
  if (is.matrix(eval(m$data,parent.frame())))  {m$data <- as.data.frame(data)}
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m,parent.frame())
  DNAME <- paste(names(mf),collapse=" by ")
  names(mf) <- NULL
  response <- attr(attr(mf,"terms"),"response")
  g <- factor(mf[[-response]])
  if (nlevels(g)!=2L) {stop("grouping factor must have exactly 2 levels")}
  DATA <- setNames(split(mf[[response]],g),c("x","y"))
  y <- do.call("fp.test",c(DATA,list(...)))
  y$data.name <- DNAME
  y
}

