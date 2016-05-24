RSS<-function(k,m,ranker){
  N<-length(ranker);
  num.samples<-m*k;
  #Each column represents one of the m*k SRS's of size k;
  SRS.index<-matrix(sample(1:N,num.samples*k),nrow=k)
  selected.rankers<-matrix(ranker[SRS.index], nrow=k)
  
  sorted<-apply(selected.rankers,2,sort)
  #Returns index of smallest, second smallest, third, and so on;
  sample.ranks<-apply(selected.rankers,2,order)
  
  output<-0
  for(i in 1:num.samples){
    index<-floor((i-1)/m)+1;
    output[i]<-SRS.index[sample.ranks[index,i],i]
  }
  
  return(sort(output))
}