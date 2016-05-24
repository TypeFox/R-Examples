#' Function to write out enrichment results
#'
#' \code{dGSEAwrite} is supposed to write out enrichment results. 
#'
#' @param eTerm an object of class "eTerm"
#' @param which_content the content will be written out. It includes two categories: i) based on "adjp" for adjusted p value, "gadjp" for globally adjusted p value, "pvalue" for p value, "FWER" for family-wise error rate, "FDR" for false discovery rate, "qvalue" for q value; ii) based on "ES" for enrichment score, "nES" for normalised enrichment score. For the former, the content is : first -1*log10-transformed, and then multiplied by -1 if nES is negative.
#' @param which_score which statistics/score will be used for declaring the significance. It can be "adjp" for adjusted p value, "gadjp" for globally adjusted p value, "FWER" for family-wise error rate, "FDR" for false discovery rate, "qvalue" for q value
#' @param cutoff a cutoff to declare the signficance. It should be used together with 'which_score'
#' @param filename a character string naming a filename
#' @param keep.significance logical to indicate whether or not to mask those insignfiicant by NA. By default, it sets to true to mask those insignfiicant by NA
#' @return 
#' a data frame with following components:
#' \itemize{
#'  \item{\code{setID}: term ID}
#'  \item{\code{setSize}: the number of genes in the set}
#'  \item{\code{name}: term name}
#'  \item{\code{namespace}: term namespace}
#'  \item{\code{distance}: term distance}
#'  \item{\code{sample names}: sample names in the next columns}
#' }
#' @note If "filename" is not NULL, a tab-delimited text file will be also written out.
#' @export
#' @seealso \code{\link{dGSEA}}
#' @include dGSEAwrite.r
#' @examples
#' #output <- dGSEAwrite(eTerm, which_content="gadjp", which_score="gadjp", filename="eTerm.txt")

dGSEAwrite <- function(eTerm, which_content=c("gadjp","adjp","pvalue","FWER","FDR","qvalue","nES","ES"), which_score=c("gadjp","adjp","FWER","FDR","qvalue","nES"), cutoff=0.1, filename=NULL, keep.significance=T) 
{

    if (class(eTerm) != "eTerm" ){
        stop("The function must apply to a 'eTerm' object.\n")
    }
    
    which_content <- match.arg(which_content)
    which_score <- match.arg(which_score)
    
    gs_size <- sapply(1:length(eTerm$gs), function(x){
        length(eTerm$gs[[x]])
    })
    
    tab <- data.frame(setID         = eTerm$set_info$setID,
                      setSize       = gs_size,
                       name         = eTerm$set_info$name,
                       namespace    = eTerm$set_info$namespace,
                       distance     = eTerm$set_info$distance
                      )
    
    ## give a sign for those scores from the statistical significance
    pos_neg <- matrix(1, nrow=nrow(eTerm$nes), ncol=ncol(eTerm$nes))
    pos_neg[eTerm$nes < 0] <- -1
    
    ## prepare content matrix
    cMatrix <- switch(which_content, ES={
        eTerm$es
    },nES={
        eTerm$nes
    }, pvalue={
        tmpD <- eTerm$pvalue
        tmpD[tmpD==0] <- 10^(-1*ceiling(max(-1*log10(tmpD[tmpD!=0])))) ## ceiling to next integer for those zeros
        -1*log10(tmpD) * pos_neg
    }, FWER={
        tmpD <- eTerm$fwer
        tmpD[tmpD==0] <- 10^(-1*ceiling(max(-1*log10(tmpD[tmpD!=0])))) ## ceiling to next integer for those zeros
        -1*log10(tmpD) * pos_neg
    }, FDR={
        tmpD <- eTerm$fdr
        tmpD[tmpD==0] <- 10^(-1*ceiling(max(-1*log10(tmpD[tmpD!=0])))) ## ceiling to next integer for those zeros
        -1*log10(tmpD) * pos_neg
    }, qvalue={
        tmpD <- eTerm$qvalue
        tmpD[tmpD==0] <- 10^(-1*ceiling(max(-1*log10(tmpD[tmpD!=0])))) ## ceiling to next integer for those zeros
        -1*log10(tmpD) * pos_neg
    }, adjp={
        tmpD <- eTerm$adjp
        tmpD[tmpD==0] <- 10^(-1*ceiling(max(-1*log10(tmpD[tmpD!=0])))) ## ceiling to next integer for those zeros
        -1*log10(tmpD) * pos_neg
    }, gadjp={
        tmpD <- eTerm$gadjp
        tmpD[tmpD==0] <- 10^(-1*ceiling(max(-1*log10(tmpD[tmpD!=0])))) ## ceiling to next integer for those zeros
        -1*log10(tmpD) * pos_neg
    }
    )
    
    ## declare significance based on which score
    flag <- switch(which_score, nES={
        which(apply(abs(eTerm$nes)>=cutoff, 1, sum) > 0)
    }, pvalue={
        tmpD <- eTerm$pvalue
        which(apply(tmpD<cutoff, 1, sum) > 0)
    }, FWER={
        tmpD <- eTerm$fwer
        which(apply(tmpD<cutoff, 1, sum) > 0)
    }, FDR={
        tmpD <- eTerm$fdr
        which(apply(tmpD<cutoff, 1, sum) > 0)
    }, qvalue={
        tmpD <- eTerm$qvalue
        which(apply(tmpD<cutoff, 1, sum) > 0)
    }, adjp={
        tmpD <- eTerm$adjp
        which(apply(tmpD<cutoff, 1, sum) > 0)
    }, gadjp={
        tmpD <- eTerm$gadjp
        which(apply(tmpD<cutoff, 1, sum) > 0)
    }
    )
    
    if(keep.significance==T){
    
        cMatrix <- switch(which_score, nES={
            cMatrix[abs(eTerm$nes)<cutoff] <- NA
            cMatrix
        }, pvalue={
            tmpD <- eTerm$pvalue
            cMatrix[tmpD>=cutoff] <- NA
            cMatrix
        }, FWER={
            tmpD <- eTerm$fwer
            cMatrix[tmpD>=cutoff] <- NA
            cMatrix
        }, FDR={
            tmpD <- eTerm$fdr
            cMatrix[tmpD>=cutoff] <- NA
            cMatrix
        }, qvalue={
            tmpD <- eTerm$qvalue
            cMatrix[tmpD>=cutoff] <- NA
            cMatrix
        }, adjp={
            tmpD <- eTerm$adjp
            cMatrix[tmpD>=cutoff] <- NA
            cMatrix
        }, gadjp={
            tmpD <- eTerm$gadjp
            cMatrix[tmpD>=cutoff] <- NA
            cMatrix
        }
        )
    }
    
    ## convert into a data frame called 'output'
    output <- cbind(tab, cMatrix)
    output <- output[flag,]
    #output <- cbind(tab[flag,], cMatrix[flag,])
    rownames(output) <- output[,1]
    
    ## If the filename is given, output data is written into a tab-delimited text file
    if(!is.null(filename)){
        utils::write.table(output, file=filename, quote=F, row.names=F, sep="\t")
        
        message(paste(c("A file called ",filename," has been successfully written!\n"), collapse=""), appendLF=T)
    }

    invisible(output)
}
