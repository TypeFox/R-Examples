################################################################
# QQ - Plot functions in package RobAStBase
################################################################



.fadeColor <- function(col,x, bg = "white"){
 ind <- seq(along=x)
 col <- .makeLenAndOrder(col,ind)
 colx <- t(sapply(ind,function(i) colorRamp(c(bg,col[i]))(x[i])))
 colv2col <- function(colvec)
   rgb(red = colvec[1], green = colvec[2], blue = colvec[3], maxColorValue = 255)
 apply(colx,1,function(x) colv2col(x))
}

## into RobAStBase
setMethod("qqplot", signature(x = "ANY",
                              y = "RobModel"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ..., distance = NormType(),
    n.adj = TRUE){

    mc <- match.call(call = sys.call(sys.parent(1)))
    if(missing(xlab)) mc$xlab <- as.character(deparse(mc$x))
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]

    if(is.null(mcl$n.CI)) mcl$n.CI <- n
    if(n.adj){
       r <- radius(neighbor(y))
       n <- floor((1-r)*n)
    }
    if(is.null(mcl$alpha.CI))
       mcl$alpha.CI <- .95
    cor <- radius(neighbor(y))
    mcl$legend.alpha <- eval(mcl$alpha.CI)
    mcl$alpha.CI <- min(eval(mcl$alpha.CI)+cor,1)


    mcl$n <- n
    mcl$y <- y@center
    mcl$legend.pref <- paste(mcl$legend.pref,"outlier-adjusted",sep="")


    xD <- fct(distance)(x)
    x.cex <- 3/(1+log(1+xD))
    mcl$cex.pch <- x.cex

    return(do.call(getMethod("qqplot", signature(x="ANY", y="ProbFamily")),
            args=mcl))
    })


## into RobAStBase
setMethod("qqplot", signature(x = "ANY", y = "InfRobModel"),
      function(x, y, n = length(x), withIdLine = TRUE, withConf = TRUE,
               withConf.pw  = withConf,  withConf.sim = withConf,
               plot.it = TRUE, xlab = deparse(substitute(x)),
               ylab = deparse(substitute(y)), ..., n.adj = TRUE){

    mc <- match.call(call = sys.call(sys.parent(1)))
    if(missing(xlab)) mc$xlab <- as.character(deparse(mc$x))
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]
    if(is.null(mcl$distance)) distance <- NormType()

    if(is.null(mcl$alpha.CI))
       mcl$alpha.CI <- .95
    cor <- radius(neighbor(y))/sqrt(n)
    mcl$legend.alpha <- eval(mcl$alpha.CI)
    mcl$alpha.CI <- min(eval(mcl$alpha.CI)+cor,1)



    if(is.null(mcl$n.CI)) mcl$n.CI <- n
    if(n.adj){
       r <- radius(neighbor(y))
       n <- floor((1-r/sqrt(n))*n)
    }
    mcl$n <- n
    mcl$y <- y@center
    mcl$legend.pref <- paste(mcl$legend.pref,"outlier-adjusted",sep="")
    
    FI <- PosSemDefSymmMatrix(FisherInfo(y@center))
    L2D <- as(diag(nrow(FI)) %*% L2deriv(y@center), "EuclRandVariable")
    L2Dx <- evalRandVar(L2D,matrix(x))[,,1]
    scx <-  solve(sqrt(FI),L2Dx)
    xD <- fct(distance)(scx)
    x.cex <- 3/(1+log(1+xD))
    mcl$cex.pch <- x.cex

    return(do.call(getMethod("qqplot", signature(x="ANY", y="ProbFamily")),
            args=mcl))
    })

## into RobAStBase
setMethod("qqplot", signature(x = "ANY",
                              y = "kStepEstimate"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...,
    exp.cex2.lbl = -.15,
    exp.cex2.pch = -.35,
    exp.fadcol.lbl = 1.85,
    exp.fadcol.pch = 1.85,
    bg = "white"
    ){

    mc <- match.call(call = sys.call(sys.parent(1)))
    if(missing(xlab)) mc$xlab <- as.character(deparse(mc$x))
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]

    IC <- pIC(y)
    if(!is(IC,"IC"))
       stop("IC of the kStepEstimator needs to be of class 'IC'")

    L2Fam <- eval(IC@CallL2Fam)

    mcl$y <- L2Fam

    if(is(IC,"HampIC")){
      dim0 <- nrow(FisherInfo(L2Fam))
      L <- as(diag(dim0)%*%L2Fam@L2deriv, "EuclRandVariable")
      L.fct <- function(x) evalRandVar(L,x)

      w.fct <- function(x)
               weight(weight(IC))(L.fct(matrix(x))[,,1])

      wx <- w.fct(x)
      mcl$order.traf <- function(x) 1/w.fct(x)

      cex.lbl <- if(is.null(mcl$cex.lbl))  par("cex")  else eval(mcl$cex.lbl)
      cex.pch <- if(is.null(mcl$cex.pch))  par("cex")  else eval(mcl$cex.pch)
      mcl$cex.lbl <- cex.lbl*wx^exp.cex2.lbl
      mcl$cex.pch <- cex.pch*wx^exp.cex2.pch

      col.lbl <- if(is.null(mcl$col.lbl))  par("col")  else eval(mcl$col.lbl)
      col.pch <- if(is.null(mcl$col.pch))  par("col")  else eval(mcl$col.pch)
      mcl$col.lbl <- .fadeColor(col.lbl,wx^exp.fadcol.lbl, bg = bg)
      mcl$col.pch <- .fadeColor(col.pch,wx^exp.fadcol.pch, bg = bg)
    }

    return(do.call(getMethod("qqplot", signature(x="ANY", y="ProbFamily")),
            args=mcl))
    })
