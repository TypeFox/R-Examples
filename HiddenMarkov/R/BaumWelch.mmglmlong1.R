BaumWelch.mmglmlong1 <- function (object, control=bwcontrol(), PSOCKcluster=NULL,
                                  tmpfile=NULL, ...){
    #   using PSOCK cluster
    tol <- control$tol
    oldLL <- -Inf
    m <- nrow(object$Pi)
    if (is.null(object$longitude)) stop("No subjects specified.")
    else {
        tmp <- table(object$longitude)
        if (min(tmp)!=max(tmp))
            stop("All subjects must have the same number of obervations.")
        N <- length(tmp)
        subnms <- names(tmp)
    }
    n <- length(object$y)/N
    #-------------------------------------------------------------
    Esteploop <- function(subnms, object, m, n){
        #  to just suppress message from mmglm1
        #  size=0 is simply to suppress an error message
        subobject <- mmglm1(NULL, object$Pi, object$delta, object$glmfamily,
                            object$beta, NULL, sigma=object$sigma,
                            nonstat=object$nonstat, size=0, msg=FALSE)
        sumcondu <- matrix(rep(0, m*n), nrow=n)
        sumcondv <- matrix(0, nrow=m, ncol=m)
        LL <- 0
        condu <- NULL
        for (subject in subnms){
            tmp <- (object$longitude==subject)
            subobject$y <- object$y[tmp]
            if (object$glmfamily$family=="binomial")
                subobject$size <- object$size[tmp]
            subobject$Xdesign <- object$Xdesign[tmp,]
            cond <- Estep.mmglm1(subobject, fortran=FALSE)
            LL <- LL + cond$LL
            sumcondu <- sumcondu + cond$u
            sumcondv <- sumcondv + apply(cond$v, MARGIN=c(2,3), FUN=sum)
            condu <- rbind(condu, cond$u)
        }
        return(list(LL=LL, sumcondu=sumcondu, sumcondv=sumcondv, condu=condu))
    }
    #-------------------------------------------------------------
    if (!is.null(PSOCKcluster)){
        numnodes <- length(PSOCKcluster)
        pernode <- trunc(N/numnodes)
        tmp <- subnms
        subnms <- list()
        for (i in 1:(numnodes-1))
            subnms[[i]] <- tmp[(1+(i-1)*pernode):(i*pernode)]
        subnms[[numnodes]] <- tmp[(1+(numnodes-1)*pernode):N]
        parallel::clusterExport(PSOCKcluster, c("mmglm1", "Estep.mmglm1", "dmmglm",
                                     "forwardback.dthmm"))
    }
    for (iter in 1:control$maxiter) {
        if (!is.null(PSOCKcluster)){
            tmp <- parallel::clusterApply(PSOCKcluster, subnms, Esteploop,
                                 object=object, m=m, n=n)
            LL <- tmp[[1]]$LL
            sumcondu <- tmp[[1]]$sumcondu
            sumcondv <- tmp[[1]]$sumcondv
            condu <- tmp[[1]]$condu
            for (i in 2:length(PSOCKcluster)){
                LL <- LL + tmp[[i]]$LL
                sumcondu <- sumcondu + tmp[[i]]$sumcondu
                sumcondv <- sumcondv + tmp[[i]]$sumcondv
                condu <- rbind(condu, tmp[[i]]$condu)
            }
        } else {
            tmp <- Esteploop(subnms, object, m, n)
            LL <- tmp$LL
            sumcondu <- tmp$sumcondu
            sumcondv <- tmp$sumcondv
            condu <- tmp$condu 
        }
        diff <- LL - oldLL
        if (control$prt) {
            cat("iter =", iter, "\n")
            cat("LL =", formatC(LL, digits=log10(1/tol)+2,
                                format="f"), "\n")
            cat("diff =", diff, "\n\n")
        }
        if (diff < 0 & control$posdiff)
            stop("Worse log-likelihood on last iteration")
        if (eval(control$converge)) break
        #----  Mstep  ----
        Pi <- diag(1/apply(sumcondv, MARGIN=1, FUN=sum)) %*% sumcondv
#  estimation of delta needs more thought
        delta <- sumcondu[1, ]/N
#       delta <- compdelta(Pi)
        tmp <- Mstep.mmglm1(object, condu)
        #-----------------
        oldLL <- LL
        object$delta <- delta
        object$Pi <- Pi
        object$beta <- tmp$beta
        object$sigma <- tmp$sigma
        if (iter %% 10){
            #   save estimates every 10th iteration
            if (!is.null(tmpfile)) save(object, file=tmpfile)
        }
    }
    rownames(object$beta) <- colnames(object$Xdesign)
    colnames(object$beta) <- paste("State", 1:length(object$delta))
    object$LL <- LL
    object$iter <- iter
    object$diff <- diff
    return(object)
}

