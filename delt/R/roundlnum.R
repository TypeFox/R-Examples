roundlnum<-function(lseq,wanted)
{
len<-length(lseq)
runner<-1
cur<-lseq[runner]
while ((cur>wanted) && (runner<len)){
    runner<-runner+1
    cur<-lseq[runner]
}
if (runner==1){
   approwanted<-cur
}
else{
  if ((wanted-cur)<=(lseq[runner-1]-wanted)){
     approwanted<-cur
  }
  else{
     approwanted<-lseq[runner-1]
  }
}

return(approwanted)
}    
