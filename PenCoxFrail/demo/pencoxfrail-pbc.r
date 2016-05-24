library(PenCoxFrail)
data("pbc")

temp <- subset(pbc, id <= 312, select=c(id:sex, stage))
pbc2 <- tmerge(temp, temp, id=id, status = event(time, status))
pbc2 <- tmerge(pbc2, pbcseq, id=id, ascites = tdc(day, ascites),
                 bili = tdc(day, bili), albumin = tdc(day, albumin),
                 protime = tdc(day, protime), alkphos = tdc(day, alk.phos))

## create two log-covariates
pbc2$log.bili <- log(pbc2$bili)
pbc2$log.protime <- log(pbc2$protime)

## cluster variable as factor
pbc2$id <- as.factor(pbc2$id)

## define death as event
pbc2$event <- as.numeric(pbc2$status==2)
  
## number of time-varying effects
m <- 2

## number of basis functions
nbasis <- 5

## number of clusters
N <- nlevels(pbc2$id)

## number of observation
n <- nrow(pbc2)

## grids for the penalty parameters
xi<-exp(seq(log(1000),log(1e-2),length=25))
zeta <- c(0,0.25,0.5,0.75,1)

## penalization of baseline hazard
penal <- 0.1

## adjust exactness of the Riemann integral approximation to the given time scale
exact <- 1

## specify all 3 formula components
fix <- as.formula("Surv(tstart, tstop, event)~1")
rnd <- list(id=~1)
vary.coef <- as.formula(" ~ log.bili + log.protime")

  y.max <- max(pbc2$time)
  W <- model.matrix(as.formula(" ~ -1+id"), pbc2)
  
  ######################
  ######## k-fold CV ###
  ######################
  set.seed(1909)
  k <- 10

  likeli <- matrix(0,length(zeta),length(xi))
  likeli.temp <- matrix(0,length(zeta),length(xi))
  
  n.cv <- rep(floor(n/k),k)
  rest <- n%%k
  
  if(rest>0){
    n.cv[1:rest] <- n.cv[1:rest] +1}
  
  which.k <- rep(1:k,n.cv)
  
  id.k <- sample(which.k,n,replace=FALSE)
  
  for(i in 1:k){
    print(paste("K.run",i,sep=": "))
    
    data.train<-pbc2[-which(id.k==i),]
    data.test<-pbc2[which(id.k==i),]
    W.test<-W[which(id.k==i),]
    
    ## calculate adaptive weights
    control <- list(start=rep(0,(m+1)*nbasis + N),q_start=1e-1,
                    smooth=list(nbasis = nbasis, penal = penal),exact=exact)
    
    cox.adapt <- try(pencoxfrail(fix=fix, vary.coef = vary.coef, rnd = rnd, 
                               data = data.train, xi = 0, control = control))
    
        
    adaptive.weights <-  cox.adapt$adaptive.weights 
    
    for(r in 1:length(zeta)){
      print(paste("zeta.run",r,sep=": "))
      
      Deltama <- as.matrix(t(rep(0,(m+1)*nbasis + N)))
      Q <- 1e-1

      for(j in 1:length(xi)){
        print(paste("xi.run",j,sep=": "))
        
        control <- list(start=Deltama[j,],q_start=Q[j],
                        smooth=list(nbasis = nbasis, penal = penal), zeta = zeta[r], 
                        exact = exact, xr = y.max)
        
        cox <- try(pencoxfrail(fix=fix, vary.coef = vary.coef, rnd = rnd, 
                            data = data.train, xi = xi[j], adaptive.weights = adaptive.weights, control = control))
        
        if(class(cox)!="try-error")
        {  
          Deltama <- rbind(Deltama,cox$Delta[cox$iter+1,])
          Q <- c(Q,cox$Q^2)

          ## prepare components necessary for calculation of the full likelihood
          del <- data.test$event; lin.pred <- W.test %*% cox$ranef
          
          U.test <- model.matrix(vary.coef, data.test)
          
          Phi.test <- bs.design(data.test$time, xl = 0, xr = y.max,
                                        spline.degree=3, nbasis = nbasis)
          Phi.test  <-cbind(Phi.test%*%cox$B.unpen.fact,Phi.test%*%cox$Pen)
          
          Phi.test.big<-matrix(0,nrow(Phi.test),(m+1)*nbasis)
          for(jj in 1:(m+1))
            Phi.test.big[,((jj-1)*nbasis+1):(jj*nbasis)] <- Phi.test * U.test[,jj]
          
          time.grid.test <- seq(0,max(data.test$time),by=cox$exact)
          
          Phi.grid <- bs.design(time.grid.test, xl = 0, xr = y.max, 
                                        spline.degree=3, nbasis = nbasis)
          Phi.grid   <-cbind(Phi.grid %*%cox$B.unpen.fact,Phi.grid %*%cox$Pen)
          
          Phi.grid.big <- t(apply(Phi.grid,1,rep,m+1))
          
          ## now calculate full likelihood (using the int.approx help function to approximate the involved integrals)
          likeli.temp[r,j] <- sum(del * (Phi.test.big %*% c(cox$baseline,cox$time.vary) + lin.pred)) - (apply(
            cbind(data.test$time,U.test),
            1,int.approx,time.grid=time.grid.test,B=Phi.grid.big,nbasis=nbasis,
            alpha=c(cox$baseline,cox$time.vary)) %*% exp(lin.pred))
        }else{
          Deltama <- rbind(Deltama,rep(NA,ncol(Deltama)))
          Q <- c(Q,NA)
          likeli.temp[r,j] <- -Inf
        }
      }}
    likeli <- likeli + likeli.temp
  }
  
  
  
  
  #############################################################################
  #############################################################################
  ### now fitting whole data set
  
  index <- which(likeli == max(likeli), arr.ind = TRUE)
  
  Deltama <- as.matrix(t(rep(0,(m+1)*nbasis + N)))
  Q <- 1e-1

  ## calculate adaptive weights
  
  control <- list(start=Deltama[1,],q_start=Q[1],
                  smooth=list(nbasis = nbasis, penal = penal),exact=exact)
  
  cox.adapt <- try(pencoxfrail(fix=fix, vary.coef = vary.coef, rnd = rnd, 
                               data = pbc2, xi = 0, control = control))
  
  for(jz in 1:index[2]){
     print(paste("xi.run",jz,sep=": "))
    
    control <- list(start=Deltama[jz,],q_start=Q[jz],
                    smooth=list(nbasis = nbasis, penal = penal), zeta =  zeta[index[1]], exact = exact)
                    
    
    cox <- try(pencoxfrail(fix=fix, vary.coef = vary.coef, rnd = rnd, 
                        data = pbc2, xi = xi[jz], adaptive.weights = cox.adapt$adaptive.weights, control = control))
    

    if(class(cox)!="try-error")
    {  
      Deltama <- rbind(Deltama,cox$Delta[cox$iter+1,])
      Q <- c(Q,cox$Q^2)
    }else{
      Deltama <- rbind(Deltama,rep(NA,ncol(Deltama)))
      Q <- c(Q,NA)
    }
  }  
  
  ###### show final fit
  
  summary(cox)
  
  #quartz()
  plot(cox)
  
  ### calculate hazard and survival for some "new" and some old subjects (with regard to id number)
  new.data <- pbc2[1:30,]
  new.data$id <- rep(c(1:2,7678:7681),table(new.data$id)[1:6])
  pred.obj <- predict(cox,newdata=new.data)
  

  ## preparations: hazard and survival can only be predicted up to the event/censoring time
  ## for this reason there are several NA's at the end of each hazard/survival curve
  notNA.ind1 <- !is.na(pred.obj$haz[,2]); notNA.ind2 <- !is.na(pred.obj$haz[,3])
  x.max <- max(c(pred.obj$time.grid[notNA.ind1],pred.obj$time.grid[notNA.ind2]))

  ## show hazard rates of an old (id=2) and a "new" (id=7678) subject
  #quartz()
  plot(pred.obj$time.grid[notNA.ind1],pred.obj$haz[notNA.ind1,2],type="l",
       ylim=c(0,0.00065),xlim=c(0,x.max),xlab="time",ylab="hazard")
  lines(pred.obj$time.grid[notNA.ind2],pred.obj$haz[notNA.ind2,3],type="l",col="red")
  legend(x=c(100,1200),y=c(0.00043,0.00055),legend=colnames(pred.obj$haz)[2:3],
         col=c("black","red"),lwd=2,lty=rep(1,2),box.lwd=2,y.intersp=1.2,cex=1.2)

  ## show survivor functions of an old (id=2) and a "new" (id=7678) subject
  #quartz()
  plot(pred.obj$time.grid[notNA.ind1],pred.obj$survival[notNA.ind1,2],type="l",
       ylim=c(0,1),xlim=c(0,x.max),xlab="time",ylab="survival")
  lines(pred.obj$time.grid[notNA.ind2],pred.obj$survival[notNA.ind2,3],type="l",col="red")
  legend(x=c(100,1200),y=c(0.4,0.6),legend=colnames(pred.obj$haz)[2:3],
         col=c("black","red"),lwd=2,lty=rep(1,2),box.lwd=2,y.intersp=1.2,cex=1.2)
  
  #############################################################################
  ######### calculate baseline hazard and all other time-varying effects ######
  #############################################################################
  
  time.seq <- seq(0,max(pbc2$tstop),by=1)

  Design<-bs.design(time.seq, xl = 0, xr = max(time.seq), spline.degree=cox$spline.degree, nbasis = cox$nbasis)
  Design  <-cbind(Design%*%cox$B.unpen.fact,Design%*%cox$Pen)
  
  smooth.effects <- Design%*%matrix(c(cox$baseline,cox$time.vary),nbasis,m+1)
  
  ## plot the baseline hazard and time-varying effects together with the results of a simple 
  ## cox model with time-constant effects using the coxph function from the survival package
  
  coxph.obj <- coxph(Surv(tstart, tstop, event) ~ log.bili + log.protime, pbc2)
  summary(coxph.obj)
  
  error <- cbind(coxph.obj$coef,sqrt(diag(coxph.obj$var)))
  
  ## calculate smooth version of baseline by the help of the gam function from the mgcv package
  detail.obj<-coxph.detail(coxph.obj)
  norm.fac<-exp(coxph.obj$coef%*%coxph.obj$mean)
  
  times<-c(0,detail.obj$t)
  hazard<-c(0,detail.obj$hazard/norm.fac)
  hazard<-c(0,hazard[-1]/diff(times)) 

  time.seq <- seq(0,max(pbc2$tstop),by=1)
  gam.data <- data.frame(time=times,haz.gam=hazard)
  names(gam.data) <- c("time","haz.gam")
  
  library(mgcv)
  gam.obj <- gam(haz.gam ~ s(time), data=gam.data)
  new.gam.data <- data.frame(time = c(time.seq,max(time.seq)+diff(time.seq)[1]))
  haz.gam.smooth<-predict(gam.obj, newdata = new.gam.data, type = "terms")[,1]+gam.obj$coef[1]
  
  
  name.vec <- c("baseline","log.bili","log.protime")
  
  #quartz()    
  par(mar=c(5, 5.5, 4, 2) + 0.1,mfrow=c(1,3))
  
  ## plot baseline hazard 
  y.lim <- c(min(c(haz.gam.smooth,exp(smooth.effects[,1]))),max(c(haz.gam.smooth,exp(smooth.effects[,1]))))
  plot(time.seq,exp(smooth.effects[,1]),type="l",ylim=y.lim,
           main=name.vec[1],xlab="time",ylab=substitute(hat(lambda)[comp](t), list(comp = 0)),
           lwd=2.5,cex.lab=1.7,cex.axis=1.5,cex.main=2) 
      lines(c(time.seq,max(time.seq)+diff(time.seq)[1]),haz.gam.smooth,col=rgb(0.8,0,0),lty=2,lwd=2.5)
      rug(jitter(pbc2$tstop[pbc2$event==1]))
      abline(h=0,lty=4)
      legend(x=c(1000,4300),y=c(y.lim[1]+0.2*y.lim[2],0.3*y.lim[2]),legend=c("coxpenfrail","coxph"),
             col=c("black",rgb(0.8,0,0)),lwd=2,lty=rep(1,2),box.lwd=2,y.intersp=1.2,cex=1.2)
      
  ## plot time-varying effects 
      for(j in 1:m)
      {
        y.lim <- c(min(c(0,error[j,1]-2*error[j,2],smooth.effects[,j+1])),max(c(error[j,1]+2*error[j,2],smooth.effects[,j+1])))
        plot(time.seq,smooth.effects[,j+1],type="l",ylim=y.lim,
           main=name.vec[j+1],xlab="time",ylab=substitute(hat(lambda)[comp](t), list(comp = j)),
           lwd=2.5,cex.lab=1.7,cex.axis=1.5,cex.main=2) 
      abline(h=error[j,1],col=rgb(0.8,0,0),lwd=2.5)
      abline(h=error[j,1]+2*error[j,2],col=rgb(0.8,0,0),lty=2,lwd=1.5)
      abline(h=error[j,1]-2*error[j,2],col=rgb(0.8,0,0),lty=2,lwd=1.5)
    abline(h=0,lty=4)
    rug(jitter(pbc2$tstop[pbc2$event==1]))
    legend(x=c(400,4300),y=c(y.lim[1]+0.15*y.lim[2],0.3*y.lim[2]),legend=c("coxpenfrail","coxph","+/- 2*SE (coxph)"),
           col=c("black",rgb(0.8,0,0),rgb(0.8,0,0)),lwd=2,lty=c(rep(1,2),2),box.lwd=2,y.intersp=1.2,cex=1.2)
    
      }
  