#' Function to obtain a list of features via splitting an input architecture
#'
#' \code{dcSplitArch} is supposed to obtain a list of features via splitting an input architecture.
#'
#' @param da an input architecture. For example, a comma-separated string
#' @param feature.mode the mode of how to define the features thereof. It can be: "supra" for combinations of one or two successive domains (including individual domains; considering the order), "individual" for individual domains only, and "comb" for all possible combinations (including individual domains; ignoring the order)
#' @param sep a character string to separate. By default, it is comma ','
#' @param ignore a character string to ignore. By default, it is '_gap_'. Ihis ignored character will affect the features defined as being 'supra' (see examples below)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a vector containing splitted features.
#' @note
#' none
#' @export
#' @seealso \code{\link{dcAlgo}}, \code{\link{dcAlgoPredict}}
#' @include dcSplitArch.r
#' @examples
#' da <- "_gap_,100895,57610,_gap_,57610,47473"
#' # get features defined as being "supra"
#' dcSplitArch(da, feature.mode="supra")
#' # get features defined as being "individual"
#' dcSplitArch(da, feature.mode="individual")
#' # get features defined as being "comb"
#' dcSplitArch(da, feature.mode="comb")

dcSplitArch <- function(da, feature.mode=c("supra","individual","comb"), sep=",", ignore="_gap_", verbose=T)
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    feature.mode <- match.arg(feature.mode)
    
    ind <- grep(ignore, da, perl=T)
    if(length(ind)!=0){
        ## first, split according to '_gap_'
        da_tmp <- unlist(strsplit(da,ignore))
        da_tmp <- da_tmp[da_tmp!='']
        pattern <- paste('^', sep, '|', sep, '$', sep='')
        da_tmp <- gsub(pattern,'', da_tmp, perl=T)
        da_tmp <- unique(da_tmp)
    }else{
        da_tmp <- da
    }
    
    if(feature.mode=='comb'){
        tmp <- paste(da_tmp, collapse=sep)
        data <- unlist(strsplit(tmp,sep))
        len <- length(data)
        
        #l <- rep(list(data), len)
        #expand.grid(l)
        
        res <- lapply(1:len, function(x){
            res <- utils::combn(data, x)
            res2 <- apply(res, 2, function(y){
                tmp <- sort(y)
                paste(tmp, collapse=sep)
            })
            unique(res2)
        })
        
    }else if(feature.mode=='supra'){
    
        res <- lapply(da_tmp, function(da){
        
            ###################
            ## then, split according to ','
            data <- unlist(strsplit(da, sep))
            len <- length(data)
    
            if(len==1){
                res <- da
                
            }else if(len <= 4){
       
                m <- sapply(0:(2^len-1),function(x){
                    as.integer(intToBits(x))
                })
                m <- t(m[1:len,])
                flag <- sapply(1:nrow(m), function(x){
                    tmp <- which(m[x,]==1)
                    if(length(tmp)==1){
                        TRUE
                    }else if(length(tmp)==0){
                        FALSE
                    }else{
                        sum(diff(tmp)==1)==(length(tmp)-1)
                    }
                })
                m <- m[flag,]
                res <- sapply(1:nrow(m), function(x){
                    paste(data[which(m[x,]==1)], collapse=',')
                })
                res <- unique(res)
        
            }else{
        
                res <- list()
                k <- 1
                st <- 1
                while(len>=st){
                    for(i in 1:(len-st+1)){
                        res[[k]] <- paste(data[i:(i+st-1)], collapse=',')
                        k <- k+1
                    }
                    st <- st+1   
                }
                res <- unique(unlist(res))
        
            }
        
            res
            ###################
        })
        
    }else if(feature.mode=='individual'){
        tmp <- paste(da_tmp, collapse=sep)
        res <- unlist(strsplit(tmp,sep))
    }
    
    res <- unique(unlist(res, use.names=F))
    
    # for SCOP fa: to remove those feature as '0'
    if(1){
        res <- res[res!='0']
    }
    
    if(verbose){
        message(sprintf("There are %d unique features (defined as being '%s').", length(res), feature.mode), appendLF=T)
    }
    
    return(res)
}
