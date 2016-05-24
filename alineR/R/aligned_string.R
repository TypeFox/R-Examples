aligned_string <-function(data,para){
  # data is a 2*n matrix containing two sequence of strings sampled from the data
  # para is a m*12 matrix contains 12 random sampled sequence. The population size of the parameter should be fixed
  string1<-vector()
  string2<-vector()
  
  ### There are mistakes in this one. After testing, it comes from the difference in the input. Maybe raw.alignment can not handle some of the strings ###
  for(i in 1:length(para[,1]))
  {
    for(j in 1:length(data[1,]))
    {
      T<-raw.alignment(c(as.character(data[1,j]),as.character(data[2,j])),para[i,1],para[i,2],para[i,3],para[i,4],para[i,5],para[i,6],para[i,7],para[i,8],para[i,9],para[i,10],para[i,11],para[i,12],para[i,13])
      string1[j]<-T$alignment1
      string2[j]<-T$alignment2
    }
    if(i==1)
    {
      sf1<-string1
      sf2<-string2
    }
    else{
      sf1<-rbind(sf1,string1)
      sf2<-rbind(sf2,string2)
    }
  }
  
 z<-list()
 z[[1]]<-sf1
 z[[2]]<-sf2
 names(z)<-c("alignment1","alignment2")
 colnames(z[[1]])<-data[1,]
 colnames(z[[2]])<-data[2,]
 return(z)
}
