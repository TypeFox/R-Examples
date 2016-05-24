

binVector <- function(piledValues, binSize) 
{
    #void binVector(int *piledValues, double *result, int *binSize, int *piledValues_Size) 
     resultSize=trunc((length(piledValues)-1)/binSize)+1

    functionReturn <- .C("C_binVector", as.integer(piledValues), res =
    double(resultSize), as.integer(binSize), as.integer(length(piledValues)))

    return(functionReturn$res)
}

#binVectorOld=function(piled, binSize)
#{
#    step=seq(1, length(piled), binSize)
#    res=vector(mode="integer", length=length(step))
#
#    index=0
#    total=length(step)
#
#    for(x in step)
#    {
#        index=index+1
#
##        cat("\n", index)
##        if((index%%trunc(total/4))==0) {cat("\n Bining",currentChr,":",round(((index)/total)*100),"%")}
##        if(index%%50000==0) {cat("\n Bining",currentChr,":",round(((index)/total)*100),"%")}
#
#        inf=x
#        sup=x+binSize-1
#
#        if(sup>length(piled))   # in case the chromosome has not a size multiple of binSize (most cases)
#                                # for the last window, we don't want to go out of bounds .
#                                # (take the end of the chrom)
#        {
#            sup=length(piled)
#        }
#
#        res[index]=mean(piled[inf:sup])
#    }
#
#    return(res)
#}
