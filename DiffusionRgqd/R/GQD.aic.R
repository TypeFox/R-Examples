GQD.aic=function(model.list,type='col')
{
  M=matrix(0,6,length(model.list))
  for(i in 1:length(model.list))
  {
    M[1,i]=model.list[[i]]$opt$convergence
    M[2,i]=length(model.list[[i]]$opt$par)
    M[3,i]=-(model.list[[i]]$opt$val)
    M[6,i]=(model.list[[i]]$model.info$N)
  }
  M[4,] = 2*M[3,]+2*M[2,]
  M[5,] = 2*M[3,]+M[2,]*log(M[6,])
  M[3:5,]=format(M[3:5,],nsmall=3)
  wh=which(as.numeric(M[4,])==min(as.numeric(M[4,])))
  M[4,wh]=paste0(' [=] ',format(as.numeric(M[4,wh]),nsmall=3))
  wh=which(as.numeric(M[5,])==min(as.numeric(M[5,])))
  M[5,wh]=paste0(' [=] ',format(as.numeric(M[5,wh]),nsmall=3))
  
  rownames(M)=c(  'Convergence      :'
                 ,'p                :'
                 ,'min likelihood)  :'
                 ,'AIC              :'
                 ,'BIC              :'
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
  res=data.frame(M)
  if(type=='col')
  {
    M=t(M)
    colnames(M)=c( 'Convergence'
                   ,'p'
                   ,'min likelihood'
                   ,'AIC'
                   ,'BIC'
                   ,'N')
    res=data.frame(M)
  }
  return(res)
}