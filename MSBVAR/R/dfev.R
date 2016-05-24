"dfev" <- function(varobj, A0=NULL, k)
{
    if (!(k>0)){
        stop("argument 'k' in dfev() must be greater than 0")
    } else if (class(varobj)==c("VAR") || class(varobj)==c("BVAR")){
        return(dfev.VAR(varobj, A0=t(chol(varobj$mean.S)), k))
    } else if (class(varobj)==c("BSVAR")){
        return(dfev.VAR(varobj, A0=solve(varobj$A0.mode), k))
    }
}

"dfev.BVAR" <- function(varobj, A0=t(chol(varobj$mean.S)), k)
{
    output <- dfev.VAR(varobj,A0,k)
    attr(output, "class") <- c("dfev")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

"dfev.BSVAR" <- function(varobj, A0=solve(varobj$A0.mode), k)
{
    output <- dfev.VAR(varobj,A0,k)
    attr(output, "class") <- c("dfev")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

"dfev.VAR" <- function(varobj, A0=t(chol(varobj$mean.S)), k)
  { m <- dim(varobj$ar.coefs)[1]
    p <- dim(varobj$ar.coefs)[3]

    # Compute the IRF
    impulses <- irf(varobj, k, A0=A0)$mhat

    # Find the cumulative innovations.  This is done by permuting the
    # array of the responses.  The irf function gives back an array
    # that m x m x k where k is the number of responses.  The last two
    # dimensions need to be permuted so we have the k x m x m arrays
    # of the responses.
    # The responses are also squared -- since they are in theory mean
    # zero, this gets them onto the right scale.
    impulses <- apply(aperm(impulses^2), c(2,3), cumsum)

    # Then compute the variances in each period.  This means we just
    # need to sum across the rows (forecast periods) for each shock
    # (now the outer array).

    var.imp <- apply(impulses,  c(1,3), sum)

    # Standardize the responses.
    for (i in 1:m)
      { impulses[,,i] <- impulses[,,i]/var.imp[,i] }

    # Scale into percentages
    errors <- 100*impulses
    # Output object
    output <- list(errors=errors, std.err=sqrt(var.imp))
    attr(output, "class") <- c("dfev")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

"print.dfev" <- function(x, latex=F, file=NULL, ...)
{   dfev.obj <- x
    errors <- dfev.obj$errors
    names <- attr(dfev.obj, "eqnames")
    std.err <- dfev.obj$std.err
    k <- dim(errors)[1]
    m <- dim(errors)[2]

    if(latex==T){
        for (i in 1:m){
            tmp <- matrix(errors[,,i], nrow=k, ncol=m)
            tmp <- cbind(std.err[,i],tmp)
            colnames(tmp) <- c("Std. Error", names)
            if(i==1){
                if(is.null(file)){
                    print(xtable(tmp, digits=rep(1,ncol(tmp)+1),
                                 caption=paste("Decomposition of Forecast Errors for a Shock to", names[i])),
                          append=F, table.placement="p")
                } else {

                      # Ensure we clobber any old file
                    print(xtable(tmp,
                                 digits=rep(1,ncol(tmp)+1),
                                 caption=paste("Decomposition of Forecast Errors for a Shock to", names[i])),
                          file=file, append=F, table.placement="p")
                }
            } else {
                if(is.null(file)){
                    print(xtable(tmp,
                                 digits=rep(1,ncol(tmp)+1),
                                 caption=paste("Decomposition of Forecast Errors for a Shock to", names[i])),
                          append=T, table.placement="p")
                } else {
                    print(xtable(tmp,
                                 digits=rep(1,ncol(tmp)+1),
                                 caption=paste("Decomposition of Forecast Errors for a Shock to", names[i])),
                          file=file, append=T, table.placement="p")
                }
            }
        }
    } else {
        for (i in 1:m){
            cat(paste("Decomposition of Forecast Errors for a Shock to", names[i], "\n"))
            cat("-------------------------------------------------------------\n")
            tmp <- matrix(errors[,,i], nrow=k, ncol=m)
            tmp <- cbind(std.err[,i],tmp)
            colnames(tmp) <- c("Std. Error", names)
            print(tmp)
            cat("-------------------------------------------------------------\n")
        }
    }
}

"summary.dfev" <- function(object, latex=F, file=NULL, ...)
{ print.dfev(object, latex=F, file=NULL, ...) }
