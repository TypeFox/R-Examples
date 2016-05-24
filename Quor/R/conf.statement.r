## File Quor/R/conf.statement.r
##
## Quor package for R (http://www.R-project.org)
## Copyright (C) 2014 Adriano Polpo, Carlos A. de B. Pereira, Cassio P. de Campos.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

conf.statement.pooled <- function(data,quantiles=NULL,ordering=NULL,verbose=TRUE,logscale=FALSE) {
    Call <- match.call()

    ## Conditions for data
    if (!is.list(data))
        stop("data is not a list")
    ngroups = length(data)
    if(verbose) cat(paste('[info] There are',ngroups,'groups\n'))

    if (is.null(quantiles)) quantiles <- c(0.5,0.5)
    if (length(quantiles) != 2)
        stop('quantiles must have two numbers for pooled confidences')
    if (is.null(ordering)) ordering <- c(1:ngroups,-(1:ngroups))
    if (!is.vector(ordering))
        stop('ordering must be a vector with the groups for pooled confidences')
    initial.time <- proc.time()[3]
    res <- c()
    for(k in ordering) {
        i <- abs(k)
        pooleddata <- NULL
        for(j in 1:ngroups) {
            if(i != j) pooleddata[[1]] <- cbind(pooleddata[[1]],data[[j]])
        }
        pooleddata[[2]] <- data[[i]]
        if(k < 0) {
            if(verbose) cat(paste('[info] CHECKING IF GROUP',i,'HAS SMALLER QUANTILE THAN POOLED GROUP OF ALL OTHERS\n'))
            cs <- conf.statement(pooleddata,quantiles,ordering=c(2,1),verbose=verbose,logscale=logscale)
        } else {
            if(verbose) cat(paste('[info] CHECKING IF GROUP',i,'HAS GREATER QUANTILE THAN POOLED GROUP OF ALL OTHERS\n'))
            cs <- conf.statement(pooleddata,quantiles,ordering=c(1,2),verbose=verbose,logscale=logscale)
        }
        res <- rbind(res,cs$confidence)
    }
    out <- NULL
    out$call <- Call
    out$pooled <- TRUE
    out$total.groups <- ngroups
    out$total.covariates <- dim(res)[2]
    out$order <- ordering
    out$quantiles <- quantiles
    out$confidence <- res
    out$run.time <- (proc.time()[3]-initial.time)
    if(verbose) cat(paste('[info] TOTAL COMPUTATION TOOK',out$run.time,'SECONDS to be completed\n'))
    class(out) <- "conf.statement"
    return(out)
}

################################################################################
## Evaluate the confidence statement for populations' quantiles
conf.statement <- function(data,quantiles=NULL,ordering=NULL,verbose=TRUE,logscale=FALSE) {
    ##    data:   is a list with all groups to be test. Each element can be
    ##            a vector with the elements, or a matrix, in which case each
    ##            row will be considered as a different covariate to be tested
    ##            All the elements in the list must have the exact same dimension
    ## quantiles: a vector of elements in [0,1] with length equal to the
    ##            number of groups. It tells us which quantile will be used
    ##            for each group. If null, then medians are compared by default.
    ## ordering:  a matrix containing one permutation of 1:n per row, where n
    ##            is the number of groups. It tells us which orderings are to be
    ##            tested. If null, then it is assumed all orderings to be tested.
    ## verbose:   if TRUE, display results on screen.
    Call <- match.call()

    ## Conditions for data
    if (!is.list(data))
        stop("data is not a list")
    ngroups = length(data)
    if(verbose) cat(paste('[info] There are',ngroups,'groups\n'))
    nrows <- 0
    if (ngroups < 2)
        stop("data must have a minimum of 2 groups")
    for (i in 1:ngroups) {
        if (is.vector(data[[i]])) {
            data[[i]] <- as.matrix(data[[i]])
            if(dim(data[[i]])[1] > 1) data[[i]] <- t(data[[i]])
        }
        if(!is.matrix(data[[i]]))
            stop(paste("data[[",i,"]] is not a vector nor matrix",sep=""))
        if (!is.numeric(data[[i]]))
            stop(paste("data[[",i,"]] is not numeric",sep=""))
        if (i==1) {
            nrows = dim(data[[i]])[1]
            if(verbose) cat(paste('[info] There are',nrows,'covariates to be analyzed\n'))
        }
        else if(nrows < 0)
            stop("Cannot mix vectors and matrices")
        if (dim(data[[i]])[2] <= 2)
            stop(paste("Sample size of data[[",i,"]] is too small (number columns <= 2)",sep=""))
        if (dim(data[[i]])[1] != nrows)
            stop("Matrices of each group must have the same number of rows")
    }
    ## Conditions for quantiles
    if (is.null(quantiles)) {
        quantiles = rep_len(0.5,ngroups) # by default, use medians
    } else {
        if (!is.vector(quantiles))
            stop("quantiles must be a vector")
        if (length(quantiles) != ngroups)
            stop("quantiles must have length equal to number of groups")
        if (length(which((quantiles <= 0) | (quantiles >= 1))) > 0)
            stop("quantiles must contain numbers strictly between 0 and 1")
    }
    ## Conditions for ordering
    if (is.null(ordering)) {
        ordering = t(matrix(data=unlist(combinat::permn(1:ngroups)),nrow=ngroups,ncol=factorial(ngroups)))
        ## by default, check all orderings
    } else {
        if (is.vector(ordering)) {
            ordering <- as.matrix(ordering)
            if(dim(ordering)[2] == 1)
                ordering <- t(ordering)
        }
        if(!is.matrix(ordering))
            stop("ordering must be a vector/matrix")
        if (dim(ordering)[2] != ngroups)
            stop("ordering must have number of rows equal to number of groups")
        ordering <- matrix(as.integer(ordering),nrow=dim(ordering)[1],ncol=dim(ordering)[2])
        ck = apply(ordering, 1, function(x) (!is.integer(x) || length(unique(x)) != ngroups || length(which((x < 1) | (x > ngroups))) > 0) )
        if(sum(ck) > 0)
            stop("ordering must contain distinct integers between 1 and number of groups")
    }
    ## Conditions to verbose
    if ((verbose != TRUE) && (verbose != FALSE))
        stop("verbose must be TRUE or FALSE")

    initial.time <- proc.time()[3]

    ## Sort data
    m <- rep_len(0,ngroups)
    for (i in 1:ngroups) {
        for(j in 1:nrows)
            data[[i]][j,] <- sort(data[[i]][j,],na.last=T)
        m[i] <- dim(data[[i]])[2]
                                        #       print(data[[i]])
        if(verbose) cat(paste('[info] Group',i,'has',m[i],'elements\n'))
    }

    ## Build cache of incomplete betas / cummulative binomials
    if(verbose) cat('[info] Building caches\n')
    ## log I_x(a,b) is pbeta(x, a, b, log.p=TRUE)
    maxm = max(m)
    lpr <- matrix(0,ngroups,maxm+1)
    upr <- matrix(-1e50,ngroups,maxm+1)
    for (i in 1:ngroups) {
        for (j in 1:m[i]) {
            ## Y_i ~ Binomial(n,p), log Pr(Y_i < j) = log Pr(Y_i <= j-1) = log F(j-1;n,p) = log I_{1-p}(n-j+1,j)
            lpr[i,j] <- pbeta(1.0-quantiles[i],m[i]-j+1,j,log.p=TRUE)
            ## Y_i ~ Binomial(n,p), log Pr(Y_i >= j) = log (1-F(j-1;n,p)) = log I_{p}(j,n-j+1)
            upr[i,j] <- pbeta(quantiles[i],j,m[i]-j+1,log.p=TRUE)
            ## cat(paste('i',i,'j',j,'lpr',lpr[i,j],'upr',upr[i,j],'\n'))
        }
    }

    ## Compute the confidence for each covariate
    if(verbose) cat('[info] Starting to process covariates\n')
    rlist <- .Call('quorccore',as.list(data),as.matrix(m),
                   as.matrix(lpr),as.matrix(upr),as.integer(ngroups),as.integer(maxm),
                   as.integer(nrows),as.matrix(ordering),as.integer(verbose),PACKAGE='Quor')
    ## break down the list with results into a matrix
    result <- t(matrix(data=unlist(rlist),nrow=nrows,ncol=dim(ordering)[1]))

####### THE FOLLOW PART HAS BEEN IMPLEMENTED IN C
    ## # R CMD SHLIB -o Quor.so core.c
    ## # dyn.load('Quor/src/Quor.so')

    ## mpr <- list()
    ## if (ngroups > 2) {
    ##     ## When more than 2 groups, build also a cache of differences in cummulative binomials
    ##     ## This takes in total O(ngroups*m^2) time, where m is the size of data for one covariate
    ##     ## This step might be slow if we think of a single covariate, but it greatly speeds up
    ##     ## computations in case there are many covariates (much more than ngroups) to be processed
    ##     mpr[[ngroups]] <- NA
    ##     for(i in 1:ngroups) {
    ##         mpr[[i]] <- matrix(-Inf,maxm+1,maxm+1)
    ##         for(j in 1:(m[i]-1)) {
    ##             for(k in (j+1):m[i]) {
    ##                 ## log( exp(lpr[i,k]) - exp(lpr(i,j)) )
    ##                 mpr[[i]][j,k] <- lpr[i,k] + log(1 - exp(lpr[i,j]-lpr[i,k]))
    ##             }
    ##         }
    ##     }
    ## }
    ## for(ires in 1:nrows) {
    ##     if(verbose && (ires %% 1000) == 0) cat('.')        
    ## ## First, for each element in a group i, mark who is the closest element
    ## ## to it in the next group i+1. The loop takes O(m) time in total. The
    ## ## allocation takes O(m*ngroups) space, which is linear in the input
    ## Jl <- matrix(NA,ngroups-1,maxm)
    ## for (i in 1:(ngroups-1)) {
    ##     jnext = 1
    ##     for(j in 1:m[i]) {
    ##         while(jnext <= m[i+1] && data[[i]][ires,j] >= data[[i+1]][ires,jnext]) jnext <- jnext+1
    ##         Jl[i,j] <- jnext
    ##     }
    ## }
    ## ## The dynamic programming begins. The allocation takes O(m*ngroups) space
    ## D <- matrix(-Inf,ngroups,maxm)
    ## ## The first line gets simply the prob of the quantile to the left of j, for each j
    ## D[1,] <- lpr[1,1:maxm]
    ## ## If there only 2 groups, next loop is skipped. If run, it takes O(m^2) time in total
    ## if(ngroups > 2) {
    ##     for (i in 2:(ngroups-1)) {
    ##         for (j in 1:m[i]) {
    ##             ## For each group, for each position j as separator, take the best solution for
    ##             ## all previous groups that are completely to the left of the separator, plus
    ##             ## the chance that the quantile of the current group will fall between that best
    ##             ## position for all previous groups and the separator j
    ##             D[i,j] <- max(D[i-1,] + mpr[[i]][Jl[i-1,],j])
    ##         }
    ##     }
    ## }
    ## ## Compose the result as the best solution found until the previous group with the prob
    ## ## that the final quantile falls to the right of the separator (the best separator is
    ## ## obtained in the maximization. The line takes O(m) time
    ## result[ires] <- max(D[ngroups-1,] + upr[Jl[ngroups-1,]])
    ##  }

    out <- NULL
    out$call <- Call
    out$pooled <- FALSE
    out$total.groups <- length(m)
    out$total.covariates <- nrows
    out$order <- ordering
    out$quantiles <- quantiles
    out$logscale <- logscale
    if(logscale)
        out$confidence <- result
    else
        out$confidence <- exp(result)
    out$run.time <- (proc.time()[3]-initial.time)
    if(verbose) cat(paste('[info] Computation took',out$run.time,'seconds to be completed\n'))
    class(out) <- "conf.statement"
    return(out)
}

print.conf.statement <- function(x, ...) {
    cat('CALL: ',paste(x$call[1]),'(',paste(x$call[-1],collapse=','),')\n')

    if (sum(x$confidence <= 1) == 0) {
        cat("It was not possible to find a confidence statement.\n")
    } else {
        cat("-----------------------------------------------------------\n")
        cat("Confidence Statement of ordered population quantiles:\n",sep="")
        cat("\n")
        if(is.matrix(x$order))
            n = dim(x$order)[1]
        else {
            if(x$pooled) n = length(x$order)
            else n = 1
        }
        cat("Number of permutations:",n,"\n")
        cat("Number of groups:",x$total.groups,"\n")
        cat("Number of variables:",x$total.covariates,"\n")
        cat("Logscale:",ifelse(x$logscale,"TRUE","FALSE"),"\n")
        cat("\n")
        dots=''
        m = dim(x$confidence)[2]
        if (m > 10) {
            dots='...'
            m = 10
        }
        conf = apply(as.matrix(x$confidence[,1:m]),2,max)
        for(i in 1:n) {
            if(x$pooled) cat(x$order[i], " : ", paste(x$confidence[i,1:m]),dots,"\n")
            else cat(paste(x$order[i,]), " : ", paste(x$confidence[i,1:m]),dots,"\n")
        }
        cat("Best : ", paste(conf),"\n")
        cat("\n")
        cat("Total time spent: ", sprintf("%.3f",x$run.time)," seconds \n",sep="")
        cat("-----------------------------------------------------------\n")
        cat("\n")
    }
}
