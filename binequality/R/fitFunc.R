fitFunc <-
function(ID, hb, bin_min, bin_max, obs_mean, ID_name, distribution=LNO,distName='LNO',links=c(muLink=identity, sigmaLink=log, nuLink=NULL, tauLink=NULL), qFunc=qLNO, quantiles=seq(0.006,0.996,length.out=1000), linksq=c(identity,exp,NULL,NULL), con=gamlss.control(c.crit = 0.1, n.cyc = 200, trace = FALSE), saveQuants=FALSE,muStart=NULL,sigmaStart=NULL, nuStart=NULL,tauStart=NULL,muFix=FALSE, sigmaFix=FALSE,nuFix=FALSE,tauFix=FALSE,freeParams=c(TRUE,TRUE,FALSE,FALSE),smartStart=FALSE, tstamp = as.numeric(Sys.time())){
  start<-Sys.time()
  
  dat<-data.frame(ID, hb, bin_min, bin_max, obs_mean)
  colnames(dat) <- c(ID_name, 'hb', 'bin_min', 'bin_max', 'obs_mean')
    
  dist<-unique(dat[,ID_name])
  
  #empty matrix to hold parameters
  parameters<-matrix(NA, ncol=4, nrow=length(dist))
  colnames(parameters)<-c('mu','sigma','nu','tau')
  rownames(parameters)<-dist
  
  if(saveQuants==TRUE){
    quants<-matrix(NA, ncol=length(quantiles), nrow=length(dist))
    colnames(quants)<-quantiles
    rownames(parameters)<-dist
  }#end if saveQuants
  
  medians<-c()
  sds<-c()
  est<-c()
  obs<-c()
  aics<-c()
  bics<-c()
  didCon<-c()
  nparams<-c()
  logLikes<-c()
  vars<-c()
  cvs<-c()
  cv2<-c()
  gini<-c()
  theil<-c()
  mld<-c()
  
  for(i in 1:length(dist)){
    if(i %% 500 == 0){print(i)}
    NAgate<-'CLOSED'
    use.i<-which(dat[,ID_name]==dist[i])
    int.i<-dat[use.i,c('bin_min','bin_max')]
    hb.i<-dat[use.i,'hb']
    intW.i<-hb.i
    rmhb<-which(intW.i==0)
    if(length(rmhb)>0){
      int.i<-int.i[-rmhb,]
      hb.i<-hb.i[-rmhb]
      intW.i<-intW.i[-rmhb]
    }#end if length(rmhb)
    if(nrow(int.i)<1){
      NAgate<-'OPEN' 
    }else{
      intCensW.i<-makeIntWeight(int.i,hb.i)
      intCens.i<-intCensW.i$int
      
      if(smartStart==TRUE& distName == "GB2"){
        mu.start.i<-as.numeric(int.i[which(hb.i==max(hb.i,na.rm=TRUE))[1],])
        mu.start.i<-mean(c(mu.start.i[1],mu.start.i[2]))
        mu.start.i<-links[[1]](mu.start.i*3.76)
        if(length(sigmaStart)==0){
          sigma.start.i<-links[[2]](2.5)
        }else{
          sigma.start.i<-sigmaStart
        }
        
        fit.i<-try(gamlss(intCens.i~1,weights=intW.i,family=cens(distribution,mu.link=links[[1]], sigma.link=links[[2]],nu.link=links[[3]], tau.link=links[[4]] ,type='interval'),mu.start=mu.start.i,sigma.start=sigma.start.i, nu.start=nuStart,tau.start=tauStart,mu.fix=muFix, sigma.fix=sigmaFix,nu.fix=nuFix,tau.fix=tauFix,control=con),silent=TRUE)
      }else{
        fit.i<-try(gamlss(intCens.i~1,weights=intW.i,family=cens(distribution,mu.link=links[[1]], sigma.link=links[[2]],nu.link=links[[3]], tau.link=links[[4]] ,type='interval'),mu.start=muStart,sigma.start=sigmaStart, nu.start=nuStart,tau.start=tauStart,mu.fix=muFix, sigma.fix=sigmaFix,nu.fix=nuFix,tau.fix=tauFix,control=con),silent=TRUE)
      }#end if/else smartStart
      
      testState<-is(fit.i)=="try-error"
      if(testState==TRUE){
        NAgate<-'OPEN' 
      }else{
        #quantiles
        quanParam.i<-getQuantilesParams(fit.i,qFunc, quantiles, linksq, freeParams,c(muStart,sigmaStart, nuStart,tauStart))
        samps.i<-quanParam.i$samps
        params.i<-quanParam.i$params
        #parameters
        parameters[i,1:length(params.i)]<-params.i
        
        #quantiles
        if(saveQuants==TRUE){
          quants[i,]<-samps.i
        }
        
        #summary stats
        medians<-c(medians,median(samps.i))
        sds<-c(sds,sd(samps.i))
        est.i<-mean(samps.i)
        var.i<-var(samps.i)
        gin.i<-giniCoef(quantiles,samps.i)
        the.i<-theilInd(samps.i)
        ll.i<-logLik(fit.i)
        est<-c(est,est.i)
        obs.i.obs<-dat[use.i[1],'hin_mean']
        mld.i<-MLD(samps.i)
      	if(length(obs.i.obs)==0){
      		obs.i.obs<-NA
      	}
      	obs<-c(obs, obs.i.obs)
        aics<-c(aics,fit.i$aic)
        bics<-c(bics,fit.i$sbc)
        didCon<-c(didCon,fit.i$converged)
        nparams<-c(nparams,attr(ll.i,'df'))
        logLikes<-c(logLikes,as.numeric(ll.i))
        vars<-c(vars,var.i)
        cvs<-c(cvs, sd(samps.i)/est.i)
        cv2<-c(cv2,(sd(samps.i)/est.i)^2)
        gini<-c(gini,gin.i)
        theil<-c(theil,the.i)
        mld<-c(mld,mld.i)
      }#end if/else testState!=TRUE
    }#end if/else nrow(int.i)<1
    if(NAgate=='OPEN'){
      medians<-c(medians,NA)
      sds<-c(sds,NA)
      ll.i<-NA
      est<-c(est,NA)
      obs.i.obs<-dat[use.i[1],'hin_mean']
      if(length(obs.i.obs)==0){
      	obs.i.obs<-NA
      }
      obs<-c(obs, obs.i.obs)
      aics<-c(aics,NA)
      bics<-c(bics,NA)
      didCon<-c(didCon,FALSE)
      nparams<-c(nparams,NA)
      logLikes<-c(logLikes,NA)
      vars<-c(vars,NA)
      cvs<-c(cvs,NA)
      cv2<-c(cv2,NA)
      gini<-c(gini,NA)
      theil<-c(theil, NA)
      mld<-c(mld,NA)
    }#end if NAgate='OPEN'
  }#end loop over dists. for i
  elasp<-Sys.time()-start
  print(elasp)
  cat('for', distName, 'fit across', length(dist), 'distributions','\n', '\n')
  distri<-rep(distName,length(est))
  datOut<-data.frame(dist,obs,distri,est,vars,cvs,cv2,gini,theil,mld,aics,bics,didCon,logLikes,nparams,medians,sds)
  colnames(datOut)<-c(ID_name,'obsMean','distribution','estMean','var','cv','cv_sqr','gini','theil','MLD','aic','bic','didConverge','logLikelihood','nparams','median','sd')
 
  
  #saving quantiles
  if(saveQuants==TRUE){
    datOut<-list('datOut' = datOut, 'timeStamp' = tstamp, 'parameters' = 'parameters', 'quantiles' = quants)
  }else{
  	datOut<-list('datOut' = datOut, 'timeStamp' = tstamp, 'parameters' = parameters, 'quantiles' = NULL)
  }
  
  return(datOut)
}
