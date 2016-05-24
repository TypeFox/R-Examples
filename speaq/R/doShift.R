doShift <-function(specSeg, shiftStep){
    nFea=length(specSeg);
    newSegment=double(nFea);
    for (j in 1:nFea)
    if (shiftStep+j>0&&shiftStep+j<=nFea) 
        newSegment[shiftStep+j]=specSeg[j];
    if (shiftStep>0){
        for (j in 1: shiftStep) newSegment[j]=newSegment[shiftStep+1];
    }
    else{
        for (j in (nFea+shiftStep): nFea) 
            newSegment[j]=newSegment[(nFea+shiftStep-1)];
    }
    return (newSegment);
}
