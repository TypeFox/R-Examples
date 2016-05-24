seqLongest <- function (
##title<< find longest sequence of TRUEs
##description<< The function determines the first position and the length of the longest
##              uninterupted series of TRUEs in a logical sequence. This for example can be
##              used to find gapless positions in a series
   array_in       ##<< logical vector: input vector
   ,na.rm=FALSE   ##<< logical: if TRUE, NAs are tretaed as FALSE
)
##\code{\link{seq.na.length}}

{
    if (na.rm)
        array_in[is.na(array_in)]=FALSE
    if (sum(array_in)==length(array_in))
    {
        pos=1
        size=length(array_in)
    } else if (sum(array_in)==0) {
        pos=0
        size=0
    } else {
         array_in_t=c(TRUE,!array_in,TRUE)
         pos    <- which(array_in_t)[which.max(diff(which(array_in_t)))]
         size   <- max(diff(which(array_in_t)))-1
    }
    #value<< list of elements
    output <- list(size=size ##<< integer: size/length of the sequence
                   ,pos=pos  ##<< integer: position/index of first element of the sequence
                   )
    output
}
