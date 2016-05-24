#Iterative quality control for daily precipitation datasets

####################
qcPrec <- function(prec,sts,inidate,enddate,parallel=TRUE,ncpu=2,printmeta=TRUE,thres=NA){

  ori=prec
  ids=colnames(prec) #station names

  #matrix of distances
  distanc=dist(cbind(sts$X,sts$Y))/1000; distanc=as.matrix(distanc)
  colnames(distanc)=sts$ID; rownames(distanc)=sts$ID

  #vector of dates
  datess=seq.Date(inidate,enddate,by='day')

  #we add a last column to insert the flag code
  prec=cbind(prec,0)

  #First function: iterative process
  qcFirst <-function(x,datess,prec,distanc,it,sts,thres){

    dw=which(x==datess)
    d=prec[dw,1:(ncol(prec)-1)]
    if(it>1 & as.numeric(prec[dw,ncol(prec)])==0) return(prec[dw,]) else {
      if(sum(is.na(d))==length(d)){
        print(paste('No data on day',x))
      } else{
        clean=data.frame(matrix(NA,ncol=4,nrow=length(d))); clean[,1] <- names(d)
        names(clean)=c('ID','obs','pred','code')
        clean$obs=as.numeric(d) #observed data in that day
        #if all data are 0...
        if(max(clean$obs,na.rm=T)==0){
          clean$pred=clean$obs
        } else{#else...

          #Now the work is by station, to that day
          for(h in 1:nrow(clean)){
            can=clean$obs[h]
            if(is.na(can)) next else{
              #neighbours
              kk=data.frame(ID=rownames(distanc),D=distanc[,which(clean$ID[h]==colnames(distanc))],
                            obs=clean$obs[match(clean$ID,rownames(distanc))],
                            stringsAsFactors=F)
              kk=kk[order(as.numeric(kk$D)),]#ordering neighbours by distance
              wna=which(is.na(kk$obs))
              if(length(wna)>0) kk=kk[-c(wna),]#and removing which have no data
              kk=kk[-c(1),]#removing first, because it is the candidate
              if(!is.na(thres)){#introduces threshold, if set
                kk=kk[-c(which(kk$D>thres)),]
              }
              if(nrow(kk)<10) {print(paste('Less than 10 nearest observations in day',
                x,'and station',clean$ID[h]))} else{ kk=kk[1:10,] #if not 10 nearest obs, go next
                  if(max(kk$obs,can)==0) next else if(max(kk$obs,na.rm=T)==0 & can>0){
                    clean$pred[h]=NA
                    clean$code[h]=1
                  } else if(min(kk$obs,na.rm=T)>0 & can==0){
                      clean$pred[h]=NA; clean$code[h]=2
                    } else{
                        #predicted
                        sts_can <- sts[which(clean$ID[h]==sts$ID),]
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

                        #more conditions...
                        if(can==0 & pb>0.5){
                          if((max((can+10)/(p+10),(p+10)/(can+10)))>10){
                            clean$pred[h]=NA
                            clean$code[h]=3
                          }
                        }
                        if(can>0){
                          if((max((can+10)/(p+10),(p+10)/(can+10)))>10){
                            clean$pred[h]=NA
                            clean$code[h]=3
                          }
                        }
                        if(can==0 & pb>0.99 & p>50){
                          clean$pred[h]=NA
                          clean$code[h]=4
                        }
                        if(can>50 & pb<0.01 & p<1){
                          clean$pred[h]=NA
                          clean$code[h]=5
                        }
                      }#next pred calculation
                    }
            }
          }#next station
        }
      }

      ww=which(!is.na(clean$code))
      if(length(ww)>0){
        k <- cbind(clean$ID[ww],rep(as.character(x),length(ww)),
                   clean$code[ww],clean$obs[ww])
        colnames(k)=c('ID','date','code','data')
        #prec[dw,unlist(lapply(k[,1],FUN=w.nam,y=colnames(prec)))] <- NA
        prec[dw,match(k[,1],colnames(prec))] <- NA
        prec[dw,ncol(prec)] <- nrow(k)
      } else prec[dw,ncol(prec)] <- 0

      return(prec[dw,])
    }
  }

  #Second function: this is like previous but using original data. Only last iteration.
  qcLast <-function(x,datess,prec,ori,distanc,sts,printmeta=printmeta,thres){
    dw=which(x==datess)
    d=prec[dw,1:(ncol(prec)-1)]
    oris=ori[dw,]#original data
    if(sum(is.na(oris))==length(d)){
      print(paste('No data on day',x))
    } else{
      clean=data.frame(matrix(NA,ncol=4,nrow=length(d))); clean[,1] <- names(d)
      names(clean)=c('ID','obs','pred','code')
      clean$obs=as.numeric(d)
      if(max(oris,na.rm=T)==0){
        clean$pred=oris
      } else{
        for(h in 1:nrow(clean)){
          can=oris[h]
          if(is.na(can)) next else{
            kk=data.frame(ID=rownames(distanc),D=distanc[,which(clean$ID[h]==colnames(distanc))],
                          obs=clean$obs[match(clean$ID,rownames(distanc))],
                          stringsAsFactors=F)
            kk=kk[order(kk$D),]
            wna=which(is.na(kk$obs))
            if(length(wna)>0) kk=kk[-c(wna),]#and removing which have no data
            kk=kk[-c(1),]
            if(!is.na(thres)){
              kk=kk[-c(which(kk$D>thres)),]
            }
            if(nrow(kk)<10) {print(paste('Less than 10 nearest observations in day',
                x,'and station',clean$ID[h]))} else {kk=kk[1:10,] #if not 10 nearest obs, go next
                  if(max(kk$obs,can)==0) next else if(max(kk$obs,na.rm=T)==0 & can>0){
                    clean$pred[h]=NA; clean$code[h]=1
                  } else if(min(kk$obs,na.rm=T)>0 & can==0){
                      clean$pred[h]=NA
                      clean$code[h]=2
                    } else{
                      sts_can <- sts[which(clean$ID[h]==sts$ID),]
                      sts_nns <- sts[match(kk$ID,sts$ID),]
                      #binomial prediction
                      b <- kk$obs; b[b>0]=1
                      DF <- data.frame(y=b,alt=sts_nns$ALT,lat=sts_nns$Y,lon=sts_nns$X)
                      fmtb <- suppressWarnings(glm(y~alt+lat+lon, data=DF, family=binomial()))
                      newdata=data.frame(alt=sts_can$ALT,lat=sts_can$Y,lon=sts_can$X)
                      pb <- predict(fmtb,newdata=newdata,type='response')
                      if(pb<0.001 & pb>0) pb <- 0.001
                      pb <- round(pb,3)
                      #data prediction
                      mini=min(kk$obs)/2
                      maxi=max(kk$obs)+(max(kk$obs)-min(kk$obs))
                      yr=as.numeric((kk$obs-mini)/(maxi-mini))
                      DF$y=yr
                      fmt <- suppressWarnings(glm(y~alt+lat+lon, data=DF, family=quasibinomial()))
                      p <- predict(fmt,newdata=newdata,type='response')
                      if(p<0.001 & p>0) p <- 0.001
                      p <- round((p*(maxi-mini))+mini,3)

                      if(can==0 & pb>0.5){
                        if((max((can+10)/(p+10),(p+10)/(can+10)))>10){
                          clean$pred[h]=NA
                          clean$code[h]=3
                        }
                      }
                      if(can>0){
                        if((max((can+10)/(p+10),(p+10)/(can+10)))>10){
                          clean$pred[h]=NA
                          clean$code[h]=3
                        }
                      }
                      if(can==0 & pb>0.99 & p>50){
                        clean$pred[h]=NA
                        clean$code[h]=4
                      }
                      if(can>50 & pb<0.01 & p<1){
                        clean$pred[h]=NA
                        clean$code[h]=5
                      }
                    }#next pred calculation
                  }
          }
        }#next station
      }
    }

    ww=which(!is.na(clean$code))
    if(length(ww)>0){
      k <- cbind(clean$ID[ww],rep(as.character(x),length(ww)),clean$code[ww],oris[ww])
      colnames(k)=c('ID','date','code','data')
      ori[dw,match(k[,1],colnames(prec))] <- NA
      if(printmeta){
        dir.create('./meta/',showWarnings = F)
        write.table(k,paste('./meta/meta_',x,'.txt',sep=''),quote=F,row.names=F,na='',sep='\t')
      }
    }
    return(ori[dw,])
  }


  #First round of QC
  if(parallel){
    sfInit(parallel=T,cpus=ncpu)
  }
  it=1
  seguir=1
  while(seguir==1){
    print(paste('Iteration',it,'of quality control'))
    #will we work with parallel code?
    if(parallel){
      d=sfLapply(datess,fun=qcFirst,datess=datess,prec=prec,distanc=distanc,it=it,sts=sts,thres=thres)
    } else {
      d=lapply(datess,FUN=qcFirst,datess=datess,prec=prec,distanc=distanc,it=it,sts=sts,thres=thres)
    }

    prec=data.frame(matrix(unlist(d),nrow=length(d),byrow=T),stringsAsFactors=F)
    colnames(prec)= ids
    gc()
    it=it+1 #go to next iteration, if needs
    if(sum(prec[,ncol(prec)])==0){
      seguir=0
    }
  }


  #Second and last round of QC
  if(parallel){
    print('Last iteration of quality control')
    d=sfLapply(datess,fun=qcLast,datess=datess,prec=prec,ori=ori,distanc=distanc,sts=sts,printmeta=printmeta,thres=thres)
  } else {
    d=lapply(datess,FUN=qcLast,datess=datess,prec=prec,ori=ori,distanc=distanc,sts=sts,printmeta=printmeta,thres=thres)
  }
  prec=data.frame(matrix(unlist(d),nrow=length(d),byrow=T),stringsAsFactors=F)
  colnames(prec)= ids
  gc()
  sfStop()
  save(prec,sts,file='cleaned.RData')
  #End of Quality Control
}
