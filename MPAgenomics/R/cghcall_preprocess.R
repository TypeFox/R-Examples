#
# This is the postsegnormalize function from CGHcall package.
#
# This function normalizes arrayCGH data after segmentation in order to find a better 0-level.
#
# @title Post-segmentation normalization
# @param segmentData  a list (see return value)
# @param inter  Interval in which the function should search for the normal level.
#
# @return a list (same format as the input list) containing
#  \describe{
#   \item{copynumber}{A matrix. Each column contains a signal of copynumber for a profile. Each row corresponds to a genomic position of a probe.}
#   \item{segmented}{A matrix of the same size as copynumber. It contains the segmented signals.}
#   \item{chromosome}{A vector of length nrow(copynumber) containing the number of the chromosome for each position.}
#   \item{startPos}{A vector of length nrow(copynumber) containing the starting genomic position of each probe.}
#   \item{featureNames}{A vector of length nrow(copynumber) containing the names of each probes.}
#   \item{sampleNames}{A vector of length ncol(copynumber) containing the names of each profiles.}
# }
# 
# @details This function recursively searches for the interval containing the most segmented data, decreasing the interval length in each recursion. 
# The recursive search makes the post-segmentation normalization robust against local maxima. 
# This function is particularly useful for profiles for which, after segmentation, the 0-level does not coincide with many segments. 
# It is more or less harmless to other profiles. We advise to keep the search interval (inter) small, in particular at the positive (gain) 
# side to avoid that the 0-level is set to a common gain level.
#
# @author Mark van de Wiel
# 
# @export
#
postsegnormalize <- function(segmentData,inter=c(-0.1,0.1))
{
    #segmentData <- seg
    seg <- segmentData$segmented
    
    values <- c()
    for (i in 1:ncol(seg)) {
                values <- c(values, median(seg[,i]));
            }    
    matrixValues    <- matrix(rep(values, nrow(seg)), ncol=ncol(seg), byrow=TRUE);
    seg <- seg - matrixValues #postseg works best when data are median normalized 
    countlevall <- apply(seg,2,function(x) {as.data.frame(table(x))})
    
    intcount <- function(int,sv){
    #int<-c(-0.5,0);sv<-segvec
        sv1 <- as.numeric(as.vector(sv[,1]))
        wh <- which(sv1<=int[2] & sv1>=int[1])
        return(sum(sv[wh,2]))
    }
    
    postsegnorm <- function(segvec,int=inter,intnr=3){
    #segvec<-countlevall[[1]];int=c(-0.30,0.1);intnr=3
        intlength <- (int[2]-int[1])/2
        gri <- intlength/intnr
        intst <- int[1]+(0:intnr)*gri
        intend <- intst+intlength
        ints <- cbind(intst,intend)
        intct <- apply(ints,1,intcount,sv=segvec)
        whmax <- which.max(intct)
        return(ints[whmax,]) 
    }
    
    postsegnorm_rec <- function(segvec,int,intnr=3){
    #segvec<-countlevall[[2]];int=c(-0.25,0.1);intnr=3
        newint <- postsegnorm(segvec,int,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        return(newint[1]+(newint[2]-newint[1])/2)
    }
    listres <- lapply(countlevall,postsegnorm_rec,int=inter)
    vecres <- c();for(i in 1:length(listres)){vecres <- c(vecres,listres[[i]])}
    
    segmentData$segmented <- t(t(seg)-vecres)
    segmentData$copynumber <- t(t(segmentData$copynumber-matrixValues)-vecres)
    return(segmentData)
}
