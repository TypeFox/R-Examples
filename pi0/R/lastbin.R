`lastbin` <-
function(p,bw=.2,trunc=TRUE){
    p=p[!is.na(p)]
    if(trunc)max(min(1,mean(p>=1-bw)/bw),0) else mean(p>=1-bw)/bw
}
