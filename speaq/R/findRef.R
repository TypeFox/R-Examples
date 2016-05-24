findRef <-function(peakList)
{
disS=matrix(data=NA,ncol=length(peakList),nrow=length(peakList));
sumDis=double(length(peakList));
for(refInd in 1:length(peakList)){
    for(tarInd in 1:length(peakList))
    if (refInd!=tarInd)
    {
        disS[refInd,tarInd]=0;
        for (i in 1:length(peakList[[tarInd]]))
            disS[refInd,tarInd]=disS[refInd,tarInd]+
                min(abs(peakList[[tarInd]][i]-peakList[[refInd]]));
    }
}

for(refInd in 1:length(peakList)){
    disS[refInd,refInd]=0;
    sumDis[refInd]=sum(disS[refInd,]);
}
orderSumdis=order(sumDis);
return(list(refInd=orderSumdis[1],orderSpec=orderSumdis));
}
