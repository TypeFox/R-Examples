#Gap filling function

gapFilling <- function(prec,sts,inidate,enddate,parallel=TRUE,ncpu=2,thres=NA){

  #matrix of distances
  distanc=dist(cbind(sts$X,sts$Y))/1000; distanc=as.matrix(distanc)
  colnames(distanc)=sts$ID; rownames(distanc)=sts$ID

  #vector of dates
  datess=seq.Date(inidate,enddate,by='day')

  fillData <-function(x,datess,prec,distanc=distanc,sts=sts,thres=thres){

    dw=which(x==datess)
    d=prec[dw,]
    if(sum(is.na(d))==length(d)){
      print(paste('No data on day',x))
      pred=data.frame(matrix(NA,ncol=7,nrow=length(d))); pred[,1] <- names(d)
      names(pred)=c('ID','obs','predb','pred1','pred2','pred3','err')
    } else{
      pred=data.frame(matrix(NA,ncol=7,nrow=length(d))); pred[,1] <- names(d)
      names(pred)=c('ID','obs','predb','pred1','pred2','pred3','err')
      pred$obs=as.numeric(d)
      if(max(pred$obs,na.rm=T)==0){
        pred$predb=0
        pred$pred1=0
        pred$pred2=0
        pred$err=0
      } else{
        for(h in 1:nrow(pred)){
          can=pred$obs[h]
          kk=data.frame(ID=rownames(distanc),D=distanc[,which(pred$ID[h]==colnames(distanc))],
                        obs=pred$obs[match(pred$ID,rownames(distanc))],
                        stringsAsFactors=F)
          kk=kk[order(kk$D),]
          wna=which(is.na(kk$obs))
          if(length(wna)>0) kk=kk[-c(wna),]#and removing which have no data
          kk=kk[-c(1),]
          if(!is.na(thres)){
            kk=kk[-c(which(kk$D>thres)),]
          }
          if(nrow(kk)<10) {kk=kk; print(paste('Less than 10 nearest observations in day',
            x,'and station',pred$ID[h]))} else {kk=kk[1:10,]
              if(max(kk$obs,na.rm=T)==0){
                pred$predb[h]=0
                pred$pred1[h]=0
                pred$pred2[h]=0
                pred$err[h]=0
              } else{
                  sts_can <- sts[which(pred$ID[h]==sts$ID),]
                  sts_nns <- sts[match(kk$ID,sts$ID),]
                  #binomial
                  b <- kk$obs; b[b>0]=1
                  DF <- data.frame(y=b,alt=sts_nns$ALT,lat=sts_nns$Y,lon=sts_nns$X)
                  fmtb <- suppressWarnings(glm(y~alt+lat+lon, data=DF, family=binomial()))
                  newdata=data.frame(alt=sts_can$ALT,lat=sts_can$Y,lon=sts_can$X)
                  pb <- predict(fmtb,newdata=newdata,type='response')
                  if(pb<0.001 & pb>0) pb <- 0.001
                  pb <- round(pb,3)

                  #data
                  mini=min(kk$obs)/2
                  maxi=max(kk$obs)+(max(kk$obs)-min(kk$obs))
                  yr=as.numeric((kk$obs-mini)/(maxi-mini))
                  DF$y=yr
                  fmt <- suppressWarnings(glm(y~alt+lat+lon, data=DF, family=quasibinomial()))
                  p <- predict(fmt,newdata=newdata,type='response')
                  if(p<0.001 & p>0) p <- 0.001
                  p <- round((p*(maxi-mini))+mini,3)

                  err <- sqrt(sum((DF$y-predict(fmt,type='response'))^2)/(length(DF$y)-3))#introduces standard error
                  if(err<0.001 & err>0) err <- 0.001
                  err <- round((err*maxi)+mini,3)
                  #asign data
                  pred$predb[h]=pb
                  pred$pred1[h]=p
                  pred$pred2[h]=p
                  if(pb<0.5) pred$pred2[h]=0
                  pred$err[h]=err
                  if(pred$predb[h]>=0.5 & pred$pred2[h]<1) pred$pred2[h]=1
                }#next pred calculation
              }
        }#next station
      }
    }
    dir.create('./days/',showWarnings = F)
    write.table(pred,paste('./days/',x,'.txt',sep=''),quote=F,row.names=F,sep='\t',na='')
  }
  #RUN gapFilling
  if(parallel){
    sfInit(parallel=T,cpus=ncpu)
  }
  if(parallel){
    print('Creating daily files')
    d=sfLapply(datess,fun=fillData,datess=datess,prec=prec,distanc=distanc,sts=sts,thres=thres)
  } else {
    d=lapply(datess,FUN=fillData,datess=datess,prec=prec,distanc=distanc,sts=sts,thres=thres)
  }
  gc()
  sfStop()

  #read predicted and standardization
  print('Re-reading data')
  aa=list.files('./days/')
  pred=matrix(NA,ncol=nrow(sts),nrow=length(datess)); colnames(pred)=sts$ID
  obs=pred
  for(i in 1:length(aa)){
    d=read.table(paste('./days/',aa[i],sep=''),header=T,sep='\t')
    obs[i,] <- d$obs
    pred[i,] <- d$pred2
  }
  print('Standardization')

  #monthly standardization
  mon1=substr(datess,1,7); mon2=substr(datess,6,7)
  for(i in 1:ncol(obs)){
    k=aggregate(obs[,i],by=list(mon1),FUN='sum',na.rm=T)
    o=aggregate(k[,2],by=list(substr(k[,1],6,7)),FUN='mean',na.rm=T)
    pp=pred[,i];pp[which(is.na(obs[,i]))]=NA
    k=aggregate(pp,by=list(mon1),FUN='sum',na.rm=T)
    p=aggregate(k[,2],by=list(substr(k[,1],6,7)),FUN='mean',na.rm=T)
    f=p[,2]/o[,2] #correction factors by months
    m=c('01','02','03','04','05','06','07','08','09','10')
    for(h in 1:12){
      w=which(m[h]==mon2)
      pred[w,i]=pred[w,i]/f[h]
    }
  }

  print('Writing final files')
  for(i in 1:length(aa)){
    d=read.table(paste('./days/',aa[i],sep=''),header=T,sep='\t')
    d$pred3 <- round(pred[i,],1) #inserting corrected values in files
    write.table(d,paste('./days/',aa[i],sep=''),quote=F,row.names=F,sep='\t',na='')
  }

  #Here we replace missing data in original by predicted and corrected values to each series
  for(i in 1:ncol(obs)){
    w=which(is.na(obs[,i]))
    if(length(w)>0) obs[w,i] <- pred[w,i]
  }

  filled=obs; rm(obs)
  save(filled,file='Filled.RData')

}
