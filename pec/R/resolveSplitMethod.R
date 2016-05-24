#' Resolve the splitMethod for estimation of prediction performance
#' 
#' The function computes a matrix of random indices obtained by drawing from
#' the row numbers of a data set either with or without replacement.  The
#' matrix can be used to repeatedly set up independent training and validation
#' sets.
#' 
#' 
#' @param splitMethod String that determines the splitMethod to use. Available
#' splitMethods are none/noPlan (no splitting), bootcv or outofbag (bootstrap
#' cross-validation), cvK (K-fold cross-validation, e.g. cv10 gives 10-fold),
#' boot632, boot632plus or boot632+, loocv (leave-one-out)
#' @param B The number of repetitions.
#' @param N The sample size
#' @param M For subsampling bootstrap the size of the subsample. Note M<N.
#' @return A list with the following components \item{name }{the official name
#' of the splitMethod} \item{internal.name }{the internal name of the
#' splitMethod} \item{index}{a matrix of indices with B columns and either N or
#' M rows, dependent on splitMethod} \item{B}{the value of the argument B}
#' \item{N}{the value of the argument N} \item{M}{the value of the argument M}
#' @author Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @keywords prediction
#' @examples
#' 
#'   # BootstrapCrossValidation: Sampling with replacement   
#'   resolvesplitMethod("BootCv",N=10,B=10)
#' 
#'   # 10-fold cross-validation: repeated 2 times
#'   resolvesplitMethod("cv10",N=10,B=2)
#' 
#'   # leave-one-out cross-validation
#'   resolvesplitMethod("loocv",N=10)
#' 
#'   resolvesplitMethod("bootcv632plus",N=10,B=2)
#'   
#'   
#' @export
resolvesplitMethod <- function(splitMethod,B,N,M){
    splitMethodName <- NULL
    k <- as.numeric(substring(grep("^cv[0-9]+$",splitMethod,value=TRUE,ignore.case=TRUE),3))
    if (length(k)==0) k <- NULL
    if (!is.null(k)){ ## classical cross-validation
        if (k==N-1){ ## this is loocv
            splitMethodName <- "LeaveOneOutCV"
            splitMethod <- "loocv"
            k <- N-1
            B <- 1
        } else{
            splitMethod <- "crossval"
            splitMethodName <- paste(k,"fold cross-validation",sep="-")
        }
    } else{
        if (length(grep("loocv",splitMethod,ignore.case=TRUE))>0){
            splitMethodName <- "LeaveOneOutCV"
            k <- N-1
            B <- 1
        }
        else{
            ## some form of bootstrap
            match.BootCv <- length(grep("boot|outofbag",splitMethod,value=FALSE,ignore.case=TRUE))>0
            if (match.BootCv==FALSE){
                splitMethod <- "noPlan"
                splitMethodName <- "no plan"
            }
            else{
                match.632 <- length(grep("632",splitMethod,value=FALSE,ignore.case=TRUE))>0
                match.plus <- length(grep("plus|\\+",splitMethod,value=FALSE,ignore.case=TRUE))>0
                if (match.632==TRUE){
                    if (match.plus==TRUE){
                        splitMethod <- "Boot632plus"
                        splitMethodName <- ".632+"
                    }
                    else{
                        splitMethod <- "Boot632"
                        splitMethodName <- ".632"
                    }
                }
                else{
                    splitMethod <- "BootCv"
                    splitMethodName <- "BootCv"}
            }
        }
    }
    if (missing(M)) M <- N
    stopifnot(M>0 && M<=N) 
    subsampling <- M!=N
    ##   if (!subsampling && resampleTraining)
    ##     stop("Resampling the training data is only available for subsampling")
    if (splitMethod %in% c("noPlan","none")) {
        B <- 0
        ##     resampleTraining <- FALSE
    }
    else{
        if (missing(B)){
            if (length(k)>0) B <- 1 # repeat k-fold CrossVal ones
            else B <- 100 # 
        }
        else if (B==0) stop("No. of resamples must be a positive integer.")
    }
    if (length(k)>0){
        if (splitMethod=="loocv")
            ResampleIndex <- data.frame(id=1:N)
        else
            ResampleIndex <- do.call("cbind",lapply(1:B,function(b){sample(rep(1:k,length.out=N))}))
    }
    else{
        if (splitMethod %in% c("Boot632plus","BootCv","Boot632")){
            ResampleIndex <- do.call("cbind",lapply(1:B,function(b){
                sort(sample(1:N,size=M,replace=!subsampling))
            }))
            colnames(ResampleIndex) <- paste("Train",1:B,sep=".")
        }
        else{
            ResampleIndex <- NULL
        }
    }
    ##   if (is.logical(resampleTraining)){
    ##     if (resampleTraining==TRUE)
    ##       resampleTrainingSize <- N
    ##   }
    ##   else{
    ##     stopifnot(resampleTraining>0 &&resampleTraining==round(resampleTraining))
    ##     resampleTrainingSize <- resampleTraining
    ##     resampleTraining <- TRUE
    ##   }
    ##   if (resampleTraining==TRUE){
    ##     if (subsampling==TRUE && resampleTrainingSize<=M)
    ##       stop("Size for resampling the training indices should exceed ",M)
    ##     ##     if (subsampling==FALSE)
    ##     ##       stop("Resampling the training indices is only allowed for subsampling")
    ##     ResampleIndex <- apply(ResampleIndex,2,function(x){
    ##       sort(c(x,sample(x,replace=TRUE,size=resampleTrainingSize-M)))
    ##     })
    ##   }
    out <- list(name=splitMethodName,
                internal.name=splitMethod,
                index=ResampleIndex,
                k=k,
                B=B,
                M=M,
                N=N)
    ##               resampleTraining=resampleTraining)
    class(out) <- "splitMethod"
    out
}

