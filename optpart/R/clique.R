clique <- function (dist,alphac,minsize=1,mult=100) 
{
     if (class(dist) != 'dist') 
         stop("The first argument must be of class 'dist'") 
     if (alphac < 0 || alphac > 1) 
         stop("alphac must be [0,1]")
     sim <- 1-as.matrix(dist)
     rows <- mult * nrow(sim)
     cols <- ncol(sim)
     top <- 0
     bottom <- 0
     orig <-0
     ds <- matrix(0,nrow=rows,ncol=cols)
     left <- rep(0,rows)
     tmp <- .Fortran('clique',
         as.double(sim),
         ds=as.integer(ds),
         as.integer(left),
         as.integer(rows),
         as.integer(cols),
         as.double(1-alphac),
         top=as.integer(top),
         bottom=as.integer(bottom),
         orig=as.integer(orig),
         PACKAGE='optpart')
    if (tmp$orig < 0) {
         print('Memory overflow.  Increase parameter mult and try again')
         out <- NULL
    } else {
         musubx <- 1 - matrix(tmp$ds,ncol=cols)[tmp$top:tmp$bottom,]
         test <- apply(musubx,1,sum) >= minsize 
         musubx <- musubx[test,]
         member <- list()
         for (i in 1:nrow(musubx)) {
           member[[i]] <- seq(1:ncol(musubx))[musubx[i,]==1]
         }
         out <- list(musubx=musubx,member=member,alphac=alphac)
         attr(out,'class') <- 'clique'
    }
    out
}

clique.test <- function (cliq,env,minsize=2,plotit=FALSE)
{
    size <- nrow(cliq$musubx)
    probs <- rep(NA,size)
    for (i in 1:size) {
        if (sum(cliq$musubx[i,]) < minsize) probs[i] <- NA
        else probs[i] <- envrtest(cliq$musubx[i,],env,plotit=plotit)$prob
        if (plotit) readline('hit return to continue')
    }
    if (plotit) {
        plot(sort(probs))
        abline(h=0.05,col=2)
    }
    invisible(probs)
}

summary.clique <- function(object,...)
{
    num <- nrow(object$musubx)
    minsize <- min(apply(object$musubx>0,1,sum))
    maxsize <- max(apply(object$musubx>0,1,sum))
    cat(paste(num,'maximal cliques at alphac = ',object$alphac),"\n")
    cat(paste('minimum size = ',minsize,"\n"))
    cat(paste('maximum size = ',maxsize,"\n"))
}

plot.clique <- function(x, panel='all', ...)
{
    if (panel == 'all' || panel == 1) {
        plot(sort(apply(x$musubx>0,1,sum)),xlab='clique',ylab='size')
        if (panel == 'all')
            readline('hit return to continue : ')
    }
    if (panel == 'all' || panel == 2) {
        plot(sort(apply(x$musubx>0,2,sum)),xlab='plot',ylab='number of cliques')
    }
}
