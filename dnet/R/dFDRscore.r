#' Function to transform fdr into scores according to log-likelihood ratio between the true positives and the false positivies and/or after controlling false discovery rate
#'
#' \code{dFDRscore} is supposed to take as input a vector of fdr, which are transformed into scores according to log-likelihood ratio between the true positives and the false positivies. Also if the FDR threshold is given, it is used to make sure that fdr below threshold are considered significant and thus scored positively. Instead, those fdr above the given threshold are considered insigificant and thus scored negatively. 
#'
#' @param fdr a vector containing a list of input fdr
#' @param fdr.threshold the given FDR threshold. By default, it is set to NULL, meaning there is no constraint. If given, those fdr with the FDR below threshold are considered significant and thus scored positively. Instead, those fdr with the FDR above given threshold are considered insigificant and thus scored negatively
#' @param scatter logical to indicate whether the scatter graph of scores against p-values should be drawn. Also indicated is the score  corresponding to the given FDR threshold (if any)
#' @return
#' \itemize{
#'  \item{\code{scores}: a vector of scores}
#' }
#' @note none
#' @export
#' @seealso \code{\link{dSVDsignif}}, \code{\link{dNetPipeline}}
#' @include dFDRscore.r
#' @examples
#' # 1) generate data with an iid matrix of 1000 x 9
#' data <- cbind(matrix(rnorm(1000*3,mean=0,sd=1), nrow=1000, ncol=3), 
#' matrix(rnorm(1000*3,mean=0.5,sd=1), nrow=1000, ncol=3), 
#' matrix(rnorm(1000*3,mean=-0.5,sd=1), nrow=1000, ncol=3))
#'
#' # 2) calculate the significance according to SVD
#' # using "fdr" significance
#' fdr <- dSVDsignif(data, signif="fdr", num.permutation=10)
#'
#' # 3) calculate the scores according to the fitted BUM and fdr=0.01
#' # no fdr threshold
#' scores <- dFDRscore(fdr)
#' # using fdr threshold of 0.01
#' scores <- dFDRscore(fdr, fdr.threshold=0.1, scatter=TRUE)

dFDRscore <- function(fdr, fdr.threshold=NULL, scatter=F)
{

    ## calculate score after taking into account the given FDR
    ## fdr below this are considered significant and thus scored positively
    ## Instead, those fdr above the given FDR are considered insigificant and thus scored negatively
    ## 
    if(is.null(fdr.threshold)){
        scores <- log2((1-fdr)/fdr)
    }else{
        if (fdr.threshold >0.5 | fdr.threshold<=0){
            warning("The function requires that the given fdr threshold falling into (0, 0.5].\n")
            return(NULL)
        }else{
            scores <- log2((1-fdr)/fdr) - log2((1-fdr.threshold)/fdr.threshold)
        }
    }
    
    ## replace those infinite values with the next finite ones
    tmp_max <- max(scores[!is.infinite(scores)])
    tmp_min <- min(scores[!is.infinite(scores)])
    scores[scores>tmp_max] <- tmp_max
    scores[scores<tmp_min] <- tmp_min

    if(scatter){
        
        grDevices::dev.new()
        
        plot(fdr, scores, pch=".", main="Scores vs FDR", xlab="FDR", ylab="Scores")
        
        if(!is.null(fdr.threshold)){
            graphics::abline(v=fdr.threshold, lty=2, col="darkgray")
            graphics::abline(h=0, lty=2, col="darkgray")
            graphics::points(fdr.threshold, 0, cex=1.5)
        }
        
    }
    
    return(scores)
}
