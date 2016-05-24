GQD.dic <-function(model.list,type='col')
{
  M=matrix(0,6,length(model.list))
  
  for(i in 1:length(model.list))
  {
    M[,i]=unlist(model.list[[i]]$model.info)[1:6]
  }
  
  whDIC = which(as.numeric(M[4,])==min(as.numeric(M[4,])))
  probs=as.numeric(M[4,whDIC])-as.numeric(M[4,])
  M[c(3,4,5),]=format(round(as.numeric(M[c(3,4,5),]),2),nsmall=3)
  wh=which(as.numeric(M[4,])==min(as.numeric(M[4,])))
  M[4,wh]=paste0(' [=] ',format(as.numeric(M[4,wh]),nsmall=3))
  
  rownames(M)=c('Elapsed time     :'
                ,'Time Homogeneous :'
                ,'p                :'
                ,'DIC              :'
                ,'pD  (effective)  :'
                ,'N                :')
  
  
  mtags=rep(0,length(model.list))
  for(i in 1:length(model.list))
  {
    mtags[i] = model.list[[i]]$model.info$Tag
  }
  
  colnames(M) = mtags
  if(all(is.na(mtags)))
  {
    colnames(M)=paste('Model',1:length(model.list))
    #warning('Some model tags are NULL!')
  }
  
  
  if(type=='row')
  {
    return(data.frame(M))
  }
  if(type=='col')
  {
    M=t(M)
    colnames(M)=c('Elapsed_Time'
                  ,'Time_Homogeneous'
                  ,'p'
                  ,'DIC'
                  ,'pD'
                  ,'N')
    
    return(data.frame(M))
  }
  
}
