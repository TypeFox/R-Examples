setMethod("getClassLabels", "Outlier", function(obj, i=1)
{
        return(which(obj@grp == levels(obj@grp)[i]))
})

setMethod("getDistance", "Outlier", function(obj)
{
   return(NULL)
})

setMethod("getFlag", "Outlier", function(obj, prob=0.975)
{
   return(obj@flag)
})

setMethod("getWeight", "Outlier", function(obj)
{
   return(obj@wt)
})

setMethod("getOutliers", "Outlier", function(obj)
{
   return(which(obj@flag == 0))
})

##
## Follow the standard methods: show, summary, plot
##
setMethod("show", "Outlier", function(object){
    cat("\nCall:\n")
    print(object@call)
    cat("-> Method: ", object@method, "\n")
    if(is.list(object@singularity))
        cat(strwrap(.MCDsingularityMsg(object@singularity, object@n.obs)), sep ="\n")

    fl <- getFlag(object)
    nout <- length(which(fl == 0))
    digits = max(3, getOption("digits") - 3)
    cat("\nNumber of outliers detected:", nout, "\n")

    print(which(fl==0))

    invisible(object)
})

##
## Follow the standard methods: show, summary, plot
##
setMethod("plot", signature(x="Outlier", y="missing"), function(x, y="missing",
                                class=1,
                                id.n=3,
                                ...){

    op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
    on.exit(par(op))

    ind <- getClassLabels(x, class)
    dist <- getDistance(x)[ind]
    const <- getCutoff(x)[class]
    flag <- getFlag(x)[ind]
    plot(dist, xlab = "Index", ylab = "Distance")
    abline(h = const)
    plot(flag, xlab = "Index", ylab = "0/1 weights", ylim = c(0, 1), ...)

})

.getGrouping <- function(grouping=NULL, n)
{
    stopifnot(!missing(n))

    if(missing(grouping) || is.null(grouping))
        grouping <- rep(0, n)

    if(length(grouping) == 1) {
        # this is the number of groups and the groups are of equal size
        ng = grouping
        ni = n/ng
        if(ng*ni < n)
            stop("nrow(x) is not divisible by the number of groups")
        grouping <- rep(0,0)
        for(i in 1:ng)
            grouping <- c(grouping, rep(i,ni))
    }else if(length(grouping) > 1 && length(grouping) < n) {
        # grouping contains a vector with the group sizes
        ng <- length(grouping)
        if(sum(grouping) != n)
            stop("nrow(x) is not equal to n1+n2+...+nn")

        gx <- rep(0,0)
        for(i in 1:ng)
            gx <- c(gx, rep(i,grouping[i]))
        grouping <- gx
    }

    if(n != length(grouping))
        stop("nrow(x) and length(grouping) are different")

    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))

    if(any(counts == 0)) {
        warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
        lev1 <- lev[counts > 0]
        g <- factor(g, levels=lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(g) <- NULL

    list(grouping=g, counts=counts, lev=lev1, ng=ng, proportions=proportions)
}
