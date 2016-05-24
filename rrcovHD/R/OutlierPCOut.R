setMethod("getCutoff", "OutlierPCOut", function(obj)
{
    const <- rep(0, length(levels(obj@grp)))
    for(i in 1:length(levels(obj@grp)))
        const[i] <- obj@covobj[[i]]$const1
    return(const)
})

setMethod("getDistance", "OutlierPCOut", function(obj)
{
    dist <- rep(0, length(obj@grp))
    for(i in 1:length(levels(obj@grp)))
    {
        class.labels <- which(obj@grp == levels(obj@grp)[i])
        dist[class.labels] <- obj@covobj[[i]]$x.dist1
    }
    return(dist)
})

##
## Follow the standard methods: show, summary, plot
##
setMethod("plot", signature(x="OutlierPCOut", y="missing"), function(x, y="missing",
                                class=1,
                                id.n=3,
                                ...){

    op <- par(mfrow=c(3,2), mar=c(4,4,2,2))
    on.exit(par(op))

    ind <- getClassLabels(x, class)

    obj.class <- x@covobj[[class]]

    n <- length(obj.class$x.dist1)
    off <- n/(n/1.5)

    # location outliers
    plot(obj.class$x.dist1, xlab="Index", ylab="Distance (location)", ...)
    if(id.n > 0)
    {
        ord <- order(obj.class$x.dist1, decreasing=TRUE)
        ii <- ord[1:id.n]
        text(x=ii+off, y=obj.class$x.dist1[ii], labels=ii)
    }
    abline(h=obj.class$const1)
    abline(h=obj.class$M1, lty=2)
    plot(obj.class$wloc, xlab="Index", ylab="Weight (location)", ylim=c(0,1), ...)
    abline(h=0)
    abline(h=1, lty=2)

    # scatter outliers
    plot(obj.class$x.dist2, xlab="Index", ylab="Distance (scatter)", ...)
    if(id.n > 0)
    {
        ord <- order(obj.class$x.dist2, decreasing=TRUE)
        ii <- ord[1:id.n]
        text(x=ii+off, y=obj.class$x.dist2[ii], labels=ii)
    }
    abline(h=obj.class$const2)
    abline(h=obj.class$M2, lty=2)
    plot(obj.class$wscat, xlab="Index", ylab="Weight (scatter)", ylim=c(0,1), ...)
    abline(h=0)
    abline(h=1,lty=2)

    # combined weights
    plot(obj.class$wfinal, xlab="Index", ylab="Weight (combined)", ylim=c(0,1), ...)
    if(id.n > 0)
    {
        ord <- order(obj.class$wfinal)
        ii <- ord[1:id.n]
        text(x=ii+off, y=obj.class$wfinal[ii], labels=ii)
    }
    abline(h=obj.class$cs)
    plot(obj.class$wfinal01, xlab="Index", ylab="Final 0/1 weight", ylim=c(0,1), ...)

})

OutlierPCOut <- function (x, ...) UseMethod("OutlierPCOut")

OutlierPCOut.formula <- function(formula, data, ..., subset, na.action)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0)
        x <- x[, -xint, drop=FALSE]
    res <- OutlierPCOut.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("OutlierPCOut")
    res@call <- cl

    res
}


OutlierPCOut.default <- function(x,
                 grouping,
                 explvar=0.99,
                 trace=FALSE,
                 ...)
{
    if(is.null(dim(x)))
        stop("x is not a matrix")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    grpx <- .getGrouping(grouping, n)
    acov <- list()
    wt <- rep(0,n)
    flag <- rep(0,n)
    outliers <- c()
    for(i in 1:grpx$ng)
    {
        class.labels <- which(grpx$grouping == grpx$lev[i])
        class <- x[class.labels,]

        xcov <- .pcout(class, explvar=explvar)

        class.outliers <- which(xcov$wfinal01 == 0)
        outliers <- c(outliers, class.labels[class.outliers])
        acov[[i]] <- xcov
        flag[class.labels] <- xcov$wfinal01
        wt[class.labels] <- xcov$wfinal
    }

    ret <- new("OutlierPCOut",
                 call=xcall,
                 counts=grpx$counts,
                 grp=grpx$grouping,
                 covobj=acov,
                 wt = wt,
                 flag=flag,
                 method="PCOUT")

    return (ret)
}

## From package mvoutlier with  minor changes
.pcout <- function(x,
                   explvar=0.99,
                   crit.M1=1/3,
                   crit.c1=2.5,
                   crit.M2=1/4,
                   crit.c2=0.99,
                   cs=0.25,
                   outbound=0.25, ...)
{

    #################################################################
    p=ncol(x)
    n=nrow(x)

    x.mad=apply(x,2,mad)
    if (any(x.mad==0))
      stop("More than 50% equal values in one or more variables!")

    #################################################################
    # PHASE 1:

    # Step 1: robustly sphere the data:
    x.sc <- scale(x,apply(x,2,median),x.mad)

    # Step 2: PC decomposition; compute p*, robustly sphere:
    ##x.wcov <- cov(x.sc)
    ##x.eig <- eigen(x.wcov)
    ##p <- ncol(x)
    ##p1 <- (1:p)[(cumsum(x.eig$val)/sum(x.eig$val)>explvar)][1]

    # or faster:
    x.svd <- svd(scale(x.sc,TRUE,FALSE))
    a <- x.svd$d^2/(n-1)
    p1 <- (1:p)[(cumsum(a)/sum(a)>explvar)][1]

    x.pc <- x.sc%*%x.svd$v[,1:p1]
    xpc.sc <- scale(x.pc,apply(x.pc,2,median),apply(x.pc,2,mad))

    # Step 3: compute robust kurtosis weights, transform to distances:
    wp <- abs(apply(xpc.sc^4,2,mean)-3)

    xpcw.sc <- xpc.sc%*%diag(wp/sum(wp))
    xpc.norm <- sqrt(apply(xpcw.sc^2,1,sum))
    x.dist1 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)

    # Step 4: determine weights according to translated biweight:
    M1 <- quantile(x.dist1,crit.M1)
    const1 <- median(x.dist1)+crit.c1*mad(x.dist1)
    w1 <- (1-((x.dist1-M1)/(const1-M1))^2)^2
    w1[x.dist1<M1] <- 1
    w1[x.dist1>const1] <- 0

    #################################################################
    # PHASE 2:

    # Step 5: compute Euclidean norms of PCs and their distances:
    xpc.norm <- sqrt(apply(xpc.sc^2,1,sum))
    x.dist2 <- xpc.norm*sqrt(qchisq(0.5,p1))/median(xpc.norm)

    # Step 6: determine weight according to translated biweight:
    M2 <- sqrt(qchisq(crit.M2,p1))
    const2 <- sqrt(qchisq(crit.c2,p1))
    w2 <- (1-((x.dist2-M2)/(const2-M2))^2)^2
    w2[x.dist2<M2] <- 1
    w2[x.dist2>const2] <- 0

    #################################################################
    # Combine PHASE1 and PHASE 2: compute final weights:
    wfinal <- (w1+cs)*(w2+cs)/((1+cs)^2)
    wfinal01 <- round(wfinal+outbound)

    list(wfinal01=wfinal01,
         wfinal=wfinal,
         wloc=w1,
         wscat=w2,
         x.dist1=x.dist1,
         x.dist2=x.dist2,
         M1=M1,
         const1=const1,
         M2=M2,
         const2=const2,
         cs=cs)
}
