#' Function to help build the contrast matrix
#'
#' \code{dContrast} is used to help build the contrast matrix
#'
#' @param level_sorted a vector of levels (usually sorted) which are contrated to each other
#' @param contrast.type the type of the contrast. It can be one of either 'average' for the contrast against the average of all levels, 'zero' for the contrast against the zero, 'sequential' for the contrast in a sequential order (it requires the levels being sorted properly), or 'pairwise' for the pairwise contrast.
#' @return 
#' a list with following components:
#' \itemize{
#'  \item{\code{each}: the contrast being specified}
#'  \item{\code{name}: the name of the contrast}
#' }
#' @note none
#' @export
#' @include dContrast.r
#' @examples
#' level_sorted <- c("L1","L2","L3","L4")
#'
#' # the contrast against the average of all levels
#' contrasts <- dContrast(level_sorted, contrast.type="average")
#'
#' # the contrast against the zero
#' contrasts <- dContrast(level_sorted, contrast.type="zero")
#'
#' # the contrast in a sequential order
#' contrasts <- dContrast(level_sorted, contrast.type="sequential")
#'
#' # the pairwise contrast
#' contrasts <- dContrast(level_sorted, contrast.type="pairwise")

dContrast <- function (level_sorted, contrast.type=c("average", "zero", "sequential", "pairwise"))
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    contrast.type <- match.arg(contrast.type)

    if(contrast.type=="average"){
        # contrast against the average
        tmp_all <- paste(level_sorted, collapse="+")
        tmp_ave <- paste("(",tmp_all,")/",length(level_sorted), sep="")
        tmp_each <- sapply(level_sorted, function(x){
            paste(x,"-",tmp_ave, sep="")
        })
        name_contrast <- sapply(level_sorted, function(x){
            paste(x, sep="")
        })
    }else if(contrast.type=="zero"){
        # contrast against the zero
        tmp_each <- sapply(level_sorted, function(x){
            paste(x, sep="")
        })
        name_contrast <- sapply(level_sorted, function(x){
            paste(x, sep="")
        })
    }else if(contrast.type=="sequential"){
        # sequential contrast
        name_contrast <- array(NA,length(level_sorted)-1)
        tmp_each <- array(NA,length(level_sorted)-1)
        for(i in 1:(length(level_sorted)-1)){
            name_contrast[i] <- paste(level_sorted[i+1],"_",level_sorted[i], sep="")
            tmp_each[i] <- paste(level_sorted[i+1],"-",level_sorted[i], sep="")
        }
    }else if(contrast.type=="pairwise"){
        # pairwise contrast
        name_contrast <- array(NA,length(level_sorted)*(length(level_sorted)-1)/2)
        tmp_each <- array(NA,length(level_sorted)*(length(level_sorted)-1)/2)
        k<-1;
        for(i in 1:(length(level_sorted)-1)){
            for(j in (i+1):length(level_sorted)){
                name_contrast[k] <- paste(level_sorted[j],"_",level_sorted[i], sep="")
                tmp_each[k] <- paste(level_sorted[j],"-",level_sorted[i], sep="")
                k <- k+1
            }
        }
    }

    contrasts <- list(each = tmp_each,
                      name = name_contrast
                 )
    
    invisible(contrasts)
}