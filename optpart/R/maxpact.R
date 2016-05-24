maxpact <- function (dist, size = NULL, alphac = NULL, mean = FALSE) 
{
    if (class(dist) != "dist") 
        stop("The first argument must be an object of class dist")
    if (is.null(size) & is.null(alphac)) stop("You must specify size or aplhac")
    if (is.null(size)) size <- attr(dist,'Size')
    sim <- 1 - as.matrix(dist)
    numplt <- nrow(sim)
    if (mean) 
        morf <- 1
    else morf <- 0
    musuby <- matrix(0, nrow = numplt, ncol = numplt)
    membry <- matrix(0, nrow = numplt, ncol = numplt)
    numset <- 0
    used <- rep(0, numplt)
    musubx <- matrix(0, nrow = numplt, ncol = numplt)
    membrx <- matrix(0, nrow = numplt, ncol = numplt)
    mnsimi <- rep(0, numplt)
    tmp <- .Fortran("maxpact", as.double(sim), as.integer(numplt), 
        as.integer(size), as.double(alphac), as.integer(morf), 
        musuby = as.double(musuby), membry = as.integer(membry), 
        numset = as.integer(numset), as.integer(used), as.double(musubx), 
        as.integer(membrx), as.double(mnsimi), PACKAGE = "optpart")
    member <- matrix(tmp$membry, nrow = numplt)
    member <- member[1:tmp$numset, 1:size]
    musuby <- matrix(tmp$musuby, nrow = numplt)
    musubx <- musuby[1:tmp$numset, 1:size]
    distname <- deparse(substitute(dist))
    if (!is.null(alphac)) {
        long <- max(apply(musubx>=alphac,1,sum))
        member[musubx<alphac] <- NA
        member <- member[,1:long]
        musubx[musubx<alphac] <- NA
        musubx <- musubx[,1:long]
    }
    out <- list(musubx = musubx, member = member, numset = tmp$numset, 
        size = size, alphac = alphac, distname = distname, 
        numele = attr(dist,'Size'))
    class(out) <- "mps"
    return(out)
}


mps.test <- function(mps,env,panel='all',main=deparse(substitute(env)),...)
{
    if (class(mps) != "mps") {
        stop("You must pass an object of class mps from maxpact()")
    }
    if (!is.numeric(env)) {
        stop("You must pass only numeric vectors as environment variables")
    }
    size <- ncol(mps$member)
    nset <- mps$numset
    null <- rep(0,nset)
    for (i in 1:nset) {
        tmp <- sample(1:length(env),size)
        nullmin <- min(env[tmp])
        nullmax <- max(env[tmp])
        null[i] <- nullmax - nullmin
    }
        
    low <- apply(mps$member,1,function(x){min(env[x])})
    high <- apply(mps$member,1,function(x){max(env[x])})
    observed <- high - low
    if (panel == 'all' | panel == 1) {
        plot(sort(null),ylim=c(0,max(null)),main=main,ylab="Within-Set Difference")
        points(sort(observed),col=2)
        if (panel == 'all')
            readline("Hit return\n")
    }
    if (panel == 'all' || panel == 2) {
        boxplot(null,observed,names=c("null","observed"),
            ylab="Within-Set Difference",main=main)
    }
    print(wilcox.test(null,observed))
}

plot.mps <- function(x, ...)
{
    plot(x$musubx[1,],ylim=c(0,1),xlab="Size",ylab="Similarity",type="n")
    for (i in 1:x$numset) {
        lines(x$musubx[i,])
    }

}

extract.mps <- function (mps) 
{
    if (class(mps) != 'mps')
        stop('You must pass an object of class mps')
    res <- list()
    for (i in 1:nrow(mps$member)) {
        tmp <- rep(FALSE,mps$numele)
        members <- c(mps$member[i,][!is.na(mps$member[i,])])
        tmp[members] <- TRUE
        res[[i]] <- tmp
    }
    res
}

