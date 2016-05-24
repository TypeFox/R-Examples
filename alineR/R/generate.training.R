# this function is creating candidate alignments for the users to rank and prepare a test set for our algorithm

generate.training<-function(raw.data,search.size=1000,table=TRUE,file.out="candidate_alignments.csv"){

  
  #raw.data should be a vector of word pairs
  raw.data<-as.matrix(raw.data)
  data_ALINE<-rbind(as.vector(encode.ALINE(raw.data[1,])),as.vector(encode.ALINE(raw.data[2,])))
  l<-length(raw.data[1,])
  
  # with the raw.data cleaned up, we move onto finding all possible alignments
  para<-random_parameter(search.size)
  align<-aligned_string(data_ALINE,para=para)
  
  # with all possible alignments generated, we now select out the identical ones
  align1<-matrix(nrow=search.size,ncol=l)
  for(i in 1:l){
    for(j in 1:(search.size)){
      align1[j,i]=paste(align[[1]][j,i],align[[2]][j,i],sep=",")
    }
  }
  
  # now we can find out all different possible alignments for different pairs
  
  alt<-list()
  for(i in 1:l){
    alt[[i]]<-list(raw.data[,i],names(summary(as.factor(align1[,i]))))
  }
  
  # now we generate the table from the result
  table.ALINE<-matrix("",nrow=300,ncol=length(raw.data[1,]))
  table.ipa<-matrix("",nrow=300,ncol=length(raw.data[1,]))
  name<-vector()
  for(i in 1:l){
    name[i]<-paste("cognate pair ",i)
    table.ipa[1:2,i]<-raw.data[,i]
    table.ALINE[1:2,i]<-data_ALINE[,i]
    for(j in 1:length(alt[[i]][[2]])){
      tem<-strsplit(alt[[i]][[2]][j],split=",")[[1]]
      table.ALINE[(2+3*j-1):(2+3*j),i]<-tem
      table.ipa[(3*j+1):(3*j+2),i]<-c(decode.ALINE(x=as.vector(raw.data[1,i]),y=tem[1]),decode.ALINE(x=as.vector(raw.data[2,i]),y=tem[2]))
     }
  }

  z.ipa<-data.frame(table.ipa)
  z.ALINE<-data.frame(table.ALINE)
  colnames(z.ipa)<-name
  colnames(z.ALINE)<-name
  
  # this is creating an output into excel files
  if(table==TRUE){
    write.csv(z.ipa,file=file.out,fileEncoding="UTF-8")
  }
  M<-list(z.ipa,z.ALINE)
  names(M)<-c("standard_ipa_symbol","ALINE_symbol")
  return(M)
  
}
