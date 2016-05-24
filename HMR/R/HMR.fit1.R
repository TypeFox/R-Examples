.HMR.fit1<-function(tid,konc,A,V,serie,ngrid,LR.always,FollowHMR,JPG,PS,PHMR,npred,xtxt,ytxt,pcttxt,MSE.zero,bracketing.tol,bracketing.maxiter,kappa.fixed)
{
  ### Stikprøvestørrelse
  n<-length(konc)
  
  ### Kammerets effektive højde
  h<-V/A

  ### HMR-analyse hvis "n>2"
  if (n>2)
  {
    ### "xOK(x)" returnerer "TRUE", hvis "x" ikke indeholder "NA", "-Inf" eller "Inf"; ellers "FALSE".
    xOK<-function(x)
    {
      if (sum(is.na(x))>0) {OK<-FALSE} else {if (max(abs(x))==Inf) {OK<-FALSE} else {OK<-TRUE}}
      OK
    }
    ### Min version af "lsfit", der afslører "collinearity", før det går galt.
    # Benytter den "qr"-test, der benyttes som standard i R - f.eks. af "lsfit"
    mylsfit<-function(x,y)
    {
      a11<-length(x); a12<-sum(x); a21<-a12; a22<-sum(x*x)
      if (xOK(c(a11,a12,a21,a22)))
      {
        XtX<-matrix(c(a11,a12,a21,a22),ncol=2,byrow=TRUE)
        dum<-try(qr(XtX)$rank,silent=TRUE)
        if ((class(dum)!='try-error')&(dum==2))
        {
          b<-as.numeric(lsfit(x,y,intercept=TRUE)$coef)
          code<-1
        } else {b<-NA; code<-0}
      } else {b<-NA; code<-0}
      list(coef=b,code=code)
    }
    ### Tester, om MSE kan beregnes.
    testMSE<-function(logkappa)
    {
      kappa<-exp(logkappa)
      x<-exp(-kappa*tid)/(-kappa*h)
      if (!xOK(x)) {ud<--1} else
      {
        dum<-mylsfit(x,konc)
        if (dum$code>0)
        {
          phi<-dum$coef[1]
          f0<-dum$coef[2]
          pkonc<-phi+f0*x
          if (!xOK(mean((konc-pkonc)^2))) {ud<--1} else {ud<-1}
        } else {ud<--1}
      }
      (ud>0)
    }
    ### MSE uden sikkerhedsnet
    MSE<-function(logkappa)
    {
      kappa<-exp(logkappa)
      x<-exp(-kappa*tid)/(-kappa*h)
      dum<-mylsfit(x,konc)$coef
      mean((konc-(dum[1]+dum[2]*x))^2)
    }
    ### MSE.list uden sikkerhedsnet
    MSE.list<-function(logkappa)
    {
      kappa<-exp(logkappa)
      x<-exp(-kappa*tid)/(-kappa*h)
      dum<-mylsfit(x,konc)$coef
      phi<-dum[1]
      f0<-dum[2]
      C0<-phi-(f0/kappa/h)
      list(C0=C0,phi=phi,f0=f0,MSE=mean((konc-(dum[1]+dum[2]*x))^2))
    }
    ### Afgrænser de mulige kappa-værdier
    indramme<-function()
    {
      # Grænser betinget af R
      logkappa.lo<-log(max(.Machine$double.xmin/h,.Machine$double.xmin))
      logkappa.up<-log(min(.Machine$double.xmax/h,.Machine$double.xmax))
      # Forsøger at gætte en "fredelig" midterværdi
      SUCCES<-TRUE
      logkappa.me<-log(1/max(tid))
      if (!testMSE(logkappa.me))
      {
        logkappa.me<-log(1/h)
        if (!testMSE(logkappa.me))
        {
          logkappa.me<-log(1)
          if (!testMSE(logkappa.me)) {SUCCES<-FALSE}
        }
      }
      # Fortsætter, hvis en midterværdi er fundet
      if (SUCCES)
      {
        # Finde "kappa.lo"
        if (!testMSE(logkappa.lo))
        {
          logk1<-logkappa.lo; logk2<-logkappa.me; logk<-(logk1+logk2)/2; iter<-1
          while ((iter<bracketing.maxiter)&(abs(logk1-logk2)>bracketing.tol)) {if (testMSE(logk)) {logk2<-logk} else {logk1<-logk}; logk<-(logk1+logk2)/2; iter<-iter+1}
          if (abs(logk1-logk2)<=bracketing.tol) {kappa.lo<-exp(logk2)} else {SUCCES<-FALSE; kappa.lo<-NA; kappa.up<-NA; code<-0}
        } else {kappa.lo<-exp(logkappa.lo)}
        # Finde "kappa.up", hvis søgningen efter "kappa.lo" lykkedes
        if (SUCCES)
        {
          if (!testMSE(logkappa.up))
          {
            logk1<-logkappa.me; logk2<-logkappa.up; logk<-(logk1+logk2)/2; iter<-1
            while ((iter<bracketing.maxiter)&(abs(logk1-logk2)>bracketing.tol)) {if (testMSE(logk)) {logk1<-logk} else {logk2<-logk}; logk<-(logk1+logk2)/2; iter<-iter+1}
            if (abs(logk1-logk2)<=bracketing.tol) {kappa.up<-exp(logk1); code<-1} else {SUCCES<-FALSE; kappa.lo<-NA; kappa.up<-NA; code<-0}
          } else {kappa.up<-exp(logkappa.up); code<-1}
        }
      } else {kappa.lo<-NA; kappa.up<-NA; code<-0}
      # Output
      list(kappa.lo=kappa.lo,kappa.up=kappa.up,code=code)
    }
    ramme<-indramme()
    # Ekstra test hvis "ramme$code=1" - findes der valide HM-modeller?
    if (ramme$code>0)
    {
      logkappa<-seq(log(ramme$kappa.lo),log(ramme$kappa.up),l=ngrid)
      vMSE<-rep(NA,ngrid); vcol<-rep(NA,ngrid)
      for (i in 1:ngrid)
      {
        dum<-MSE.list(logkappa[i])
        vMSE[i]<-dum$MSE
        if ((dum$C0<=0)|(dum$phi<=0)) {vcol[i]<-2} else {vcol[i]<-3}
      }
      if (sum(vcol==3)==0) {ramme$code<-0}
    }
  } else {ramme<-list(code=0)}

  ### Valg af dataanalyse
  # Hvis HMR kan anvendes
  if (ramme$code>0)
  {
    # Global gittersøgning
    logkappa<-seq(log(ramme$kappa.lo),log(ramme$kappa.up),l=ngrid)
    vMSE<-rep(NA,ngrid); vcol<-rep(NA,ngrid)
    for (i in 1:ngrid)
    {
      dum<-MSE.list(logkappa[i])
      vMSE[i]<-dum$MSE
      if ((dum$C0<=0)|(dum$phi<=0)) {vcol[i]<-2} else {vcol[i]<-3}
    }
    # Betinget gittersøgning
    dum<-1:length(vMSE); ok<-dum[vcol==3]; bmin.vMSE<-min(vMSE[ok])
    kand<-dum[ok][vMSE[ok]<=bmin.vMSE]
    if (length(kand)>1)
    {
      # Best-in-grid estimatet er ikke entydigt bestemt; vælger største (~No flux) eller mindste (~LR) kandidat
      if (abs(vMSE[ok][1]-bmin.vMSE)<abs(vMSE[ok][length(vMSE[ok])]-bmin.vMSE)) {big<-kand[1]} else {big<-kand[length(kand)]}
    } else {big<-kand}
    logkappa.big<-logkappa[big]
    MSE.big<-MSE(logkappa.big)
    # Betinget optimering
    if ((big==1)|(big==ngrid))
    {
      logkappa.opt<-logkappa.big
      MSE.opt<-MSE.big
    } else
    {
      # Find interval omkring "logkappa[big]" indenfor gitteret
      lo.i<-big
      GREEN<-TRUE
      while ((lo.i>1)&(GREEN)) {lo.i<-lo.i-1; GREEN<-(vcol[lo.i]==3)}
      lo<-logkappa[lo.i+1]
      up.i<-big
      GREEN<-TRUE
      while ((up.i<ngrid)&(GREEN)) {up.i<-up.i+1; GREEN<-(vcol[up.i]==3)}
      up<-logkappa[up.i-1]
      # Optimere
      opt<-optimize(f=MSE,lower=lo,upper=up)
      logkappa.opt<-opt$minimum; MSE.opt<-opt$objective
      if (MSE.big<MSE.opt)
      {
        opt<-optimize(f=MSE,lower=max(logkappa.big-log(2),lo),upper=min(logkappa.big+log(2),up))
        logkappa.opt<-opt$minimum; MSE.opt<-opt$objective
        if (MSE.big<MSE.opt) {logkappa.opt<-logkappa.big; MSE.opt<-MSE.big}
      }
      # Kontrollere at optimum svarer til valid HM-model
      EST<-MSE.list(logkappa.opt)  
      if ((EST$C0<=0)|(EST$phi<=0)) {logkappa.opt<-logkappa.big; MSE.opt<-MSE.big}
    }
    # Test for linearitet og manglende fit
    maintext<-'HMR' # Det er kontrolleret ovenfor, at "logkappa.opt" svarer til valid HM-model
    if (((vMSE[vcol==3][1]-MSE.big)<=MSE.zero)&((vMSE[vcol==3][length(vMSE[vcol==3])]-MSE.big)>MSE.zero))
    {
      EST<-as.numeric(lsfit(tid,konc)$coef)
      if (EST[1]>0) {maintext<-'LR'} else {maintext<-'HMR'}
    }
    if ((vMSE[vcol==3][length(vMSE[vcol==3])]-MSE.big)<=MSE.zero) {maintext<-'No flux/LR'}
    # Valg af analyse
    if (FollowHMR)
    { # Bare kør løs!
      if (maintext=='HMR') {dacode<-1} else
        if (maintext=='LR') {dacode<-2} else {dacode<-3}
    } else
    { # Valg af analyse foretages af brugeren
      options(warn=-1)
      par(mfrow=c(2,2),oma=c(0,0,2,0),bty='n',pty='m')
      plot(logkappa,vMSE,type='p',pch=16,col=vcol,xlab=expression(paste('Log(',kappa,')',sep='')),
      ylab='MSE',main='MSE criterion')
      points(logkappa.opt,MSE.opt,type='p',pch=15,cex=1.5,col=4)
      lo<-max(logkappa.opt-log(2),log(ramme$kappa.lo))
      up<-min(logkappa.opt+log(2),log(ramme$kappa.up))
      x<-logkappa[(logkappa>=lo)&(logkappa<=up)]
      if (length(x)>0)
      {
        y<-vMSE[(logkappa>=lo)&(logkappa<=up)]
        plot(x,y,type='p',pch=16,col=vcol[(logkappa>=lo)&(logkappa<=up)],xlab=expression(paste('Log(',kappa,')',sep='')),
        ylab='MSE',main='MSE criterion - zoom on HMR optimum')
        points(logkappa.opt,MSE.opt,type='p',pch=15,cex=1.5,col=4)
      } else
      {
        plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='',main='MSE criterion - zoom on HMR optimum')
        lines(c(0,1),c(0,0),lty=1,lwd=1,col=1)
        lines(c(0,1),c(1,1),lty=1,lwd=1,col=1)
        lines(c(0,0),c(0,1),lty=1,lwd=1,col=1)
        lines(c(1,1),c(0,1),lty=1,lwd=1,col=1)
        text(0.5,0.5,'Too few grid-points!',adj=0.5,col=2)
      }
      ptid<-seq(0,max(tid),l=npred)
      kappa<-exp(logkappa.opt); EST<-MSE.list(logkappa.opt); phi<-EST$phi; f0<-EST$f0; C0<-EST$C0
      x<-exp(-kappa*ptid)/(-kappa*h); pkonc.HMR<-phi+f0*x
      maxy<-min(max(konc,C0),2*max(konc))
      plot(c(0,max(tid)),c(min(konc,C0),maxy),type='n',cex=1.5,xlab=xtxt,ylab=ytxt,main=paste('Recommendation: ',maintext,sep=''))
      points(tid,konc,type='p',pch=16,cex=1.5)
      EST<-as.numeric(lsfit(tid,konc)$coef)
      pkonc.LR<-EST[1]+EST[2]*ptid
      pkonc.NE<-rep(mean(konc),length(ptid))
      lines(ptid,pkonc.HMR,lty=1,col=4)
      lines(ptid,pkonc.LR,lty=1,col=colors()[498])
      lines(ptid,pkonc.NE,lty=2,col=1)
      if (C0>max(konc)) {legend(median(tid),maxy,legend=c('HMR','LR','No flux'),lty=c(1,1,2),col=c(4,colors()[498],1))} else
      {legend(0,max(konc),legend=c('HMR','LR','No flux'),lty=c(1,1,2),col=c(4,colors()[498],1))}
      plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='',main='')
      polygon(c(0.1,0.1,0.45,0.45,0.1),c(0.7,0.9,0.9,0.7,0.7),density=-1,col=4,border=1,lty=1,lwd=1)
      text(0.275,0.8,'HMR',adj=0.5,col=colors()[1],cex=1.5)
      polygon(c(0.1,0.1,0.45,0.45,0.1),c(0.4,0.6,0.6,0.4,0.4),density=-1,col=colors()[498],border=1,lty=1,lwd=1)
      text(0.275,0.5,'LR',adj=0.5,col=1,cex=1.5)
      polygon(c(0.1,0.1,0.45,0.45,0.1),c(0.1,0.3,0.3,0.1,0.1),density=-1,col=1,border=1,lty=1,lwd=1)
      text(0.275,0.2,'No flux',adj=0.5,col=colors()[1],cex=1.5)
      polygon(c(0.55,0.55,0.9,0.9,0.55),c(0.1,0.3,0.3,0.1,0.1),density=-1,col=2,border=1,lty=1,lwd=1)
      text(0.725,0.2,'CANCEL',adj=0.5,col=colors()[1],cex=1.5)
      mtext(paste('Data series: ',serie,pcttxt,sep=''),side=3,line=0,outer=TRUE,cex=1.5)
      VALGT<-FALSE
      while (!VALGT)
      {
        pkt<-locator(1)
        if ((pkt[1]>=0.1)&(pkt[1]<=0.45)&(pkt[2]>=0.7)&(pkt[2]<=0.9)) {dacode<-1; VALGT<-TRUE} # HMR
        if ((pkt[1]>=0.1)&(pkt[1]<=0.45)&(pkt[2]>=0.4)&(pkt[2]<=0.6)) {dacode<-2; VALGT<-TRUE} # LR
        if ((pkt[1]>=0.1)&(pkt[1]<=0.45)&(pkt[2]>=0.1)&(pkt[2]<=0.3)) {dacode<-3; VALGT<-TRUE} # No flux
        if ((pkt[1]>=0.55)&(pkt[1]<=0.9)&(pkt[2]>=0.1)&(pkt[2]<=0.3)) {dacode<-0; VALGT<-TRUE} # CANCEL
      }
      options(warn=0)
      # JPG, PS
      dum<-c(JPG,PS)
      for (fig in 1:2)
        if (dum[fig])
        {
          if (fig==1)
            jpeg(filename=paste(serie,'.jpg',sep=''),width=800,height=800,pointsize=12,quality=100)
          else
            postscript(file=paste(serie,'.ps',sep=''),horizontal=TRUE)
          par(mfrow=c(1,1),oma=c(0,0,0,0),bty='n',pty='m')
          maxy<-min(max(konc,C0),2*max(konc))
          plot(c(0,max(tid)),c(min(konc,C0),maxy),type='n',cex=1.5,xlab=xtxt,ylab=ytxt,main=paste(serie,': ',maintext,sep=''))
          points(tid,konc,type='p',pch=16,cex=1.5)
          lines(ptid,pkonc.HMR,lty=1,col=4)
          lines(ptid,pkonc.LR,lty=1,col=colors()[498])
          lines(ptid,pkonc.NE,lty=2,col=1)
          if (C0>max(konc)) {legend(median(tid),maxy,legend=c('HMR','LR','No flux'),lty=c(1,1,2),col=c(4,colors()[498],1))} else
          {legend(0,max(konc),legend=c('HMR','LR','No flux'),lty=c(1,1,2),col=c(4,colors()[498],1))}
          dev.off()
          par(mfrow=c(2,2),oma=c(0,0,2,0),bty='n',pty='m')
        }
    }
  } else
  # Hvis HMR ikke kan anvendes
  {
    # Valg af analyse
    if (FollowHMR)
    { # Forsigtighedsprincip: No flux (dacode=3)
      dacode<-3
    } else
    { # Valg af analyse foretages af brugeren
      pframe<-function() {lines(c(0,1),c(0,0),lty=1,lwd=2,col=1); lines(c(0,1),c(1,1),lty=1,lwd=2,col=1); lines(c(0,0),c(0,1),lty=1,lwd=2,col=1); lines(c(1,1),c(0,1),lty=1,lwd=2,col=1)}
      par(mfrow=c(2,2),oma=c(0,0,2,0),bty='n',pty='m')
      plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='',main='MSE criterion'); pframe()
      plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='',main='MSE criterion - zoom on HMR optimum'); pframe()
      plot(c(0,max(tid)),c(min(konc),max(konc)),type='n',cex=1.5,xlab=xtxt,ylab=ytxt,main='HMR could not be applied')
      points(tid,konc,type='p',pch=16,cex=1.5)
      ptid<-seq(0,max(tid),l=npred)
      EST<-as.numeric(lsfit(tid,konc)$coef)
      pkonc.LR<-EST[1]+EST[2]*ptid
      pkonc.NE<-rep(mean(konc),length(ptid))
      lines(ptid,pkonc.LR,lty=1,col=colors()[498])
      lines(ptid,pkonc.NE,lty=2,col=1)
      legend(0,max(konc),legend=c('LR','No flux'),lty=c(1,2),col=c(colors()[498],1))
      plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='',main='')
      polygon(c(0.1,0.1,0.45,0.45,0.1),c(0.4,0.6,0.6,0.4,0.4),density=-1,col=colors()[498],border=1,lty=1,lwd=1)
      text(0.275,0.5,'LR',adj=0.5,col=1,cex=1.5)
      polygon(c(0.1,0.1,0.45,0.45,0.1),c(0.1,0.3,0.3,0.1,0.1),density=-1,col=1,border=1,lty=1,lwd=1)
      text(0.275,0.2,'No flux',adj=0.5,col=colors()[1],cex=1.5)
      polygon(c(0.55,0.55,0.9,0.9,0.55),c(0.1,0.3,0.3,0.1,0.1),density=-1,col=2,border=1,lty=1,lwd=1)
      text(0.725,0.2,'CANCEL',adj=0.5,col=colors()[1],cex=1.5)
      mtext(paste('Data series: ',serie,pcttxt,sep=''),side=3,line=0,outer=TRUE,cex=1.5)
      VALGT<-FALSE
      while (!VALGT)
      {
        pkt<-locator(1)
        if ((pkt[1]>=0.1)&(pkt[1]<=0.45)&(pkt[2]>=0.4)&(pkt[2]<=0.6)) {dacode<-2; VALGT<-TRUE} # LR
        if ((pkt[1]>=0.1)&(pkt[1]<=0.45)&(pkt[2]>=0.1)&(pkt[2]<=0.3)) {dacode<-3; VALGT<-TRUE} # No flux
        if ((pkt[1]>=0.55)&(pkt[1]<=0.9)&(pkt[2]>=0.1)&(pkt[2]<=0.3)) {dacode<-0; VALGT<-TRUE} # CANCEL
      }
      # JPG, PS
      dum<-c(JPG,PS)
      for (fig in 1:2)
        if (dum[fig])
        {
          if (fig==1)
            jpeg(filename=paste(serie,'.jpg',sep=''),width=800,height=800,pointsize=12,quality=100)
          else
            postscript(file=paste(serie,'.ps',sep=''),horizontal=TRUE)
          par(mfrow=c(1,1),oma=c(0,0,0,0),bty='n',pty='m')
          plot(c(0,max(tid)),c(min(konc),max(konc)),type='n',cex=1.5,xlab=xtxt,ylab=ytxt,main=paste(serie,': HMR could not be applied',sep=''))
          points(tid,konc,type='p',pch=16,cex=1.5)
          lines(ptid,pkonc.LR,lty=1,col=colors()[498])
          lines(ptid,pkonc.NE,lty=2,col=1)
          legend(0,max(konc),legend=c('LR','No flux'),lty=c(1,2),col=c(colors()[498],1))
          dev.off()
          par(mfrow=c(2,2),oma=c(0,0,2,0),bty='n',pty='m')
        }
    }
  }

  ### Datanalyse
  PHMR.ptid<-NA; PHMR.pkonc<-NA
  if (dacode>0)
  {
    if (dacode==1)
    {
      # HMR
      kappa<-exp(logkappa.opt)
      x<-exp(-kappa*tid)/(-kappa*h)
      dum<-lsfit(x,konc)
      phi<-as.numeric(dum$coef)[1]
      f0.est<-as.numeric(dum$coef)[2]
      if (PHMR) {PHMR.ptid<-seq(0,max(tid),l=npred); PHMR.pkonc<-phi+f0.est*(exp(-kappa*PHMR.ptid)/(-kappa*h))}
      if (kappa.fixed)
      {
        if (n>3)
        {
          f0.se<-as.numeric(ls.diag(dum)$std.err)[2]
          f0.p<-2*pt(q=-abs(f0.est/f0.se),df=n-2)
          fraktil<-qt(p=0.975,df=n-2)
          f0.lo95<-f0.est-fraktil*f0.se
          f0.up95<-f0.est+fraktil*f0.se
          method<-'HMR'
          advarsel<-'None'
        } else
        {
          f0.p<-NA; f0.lo95<-NA; f0.up95<-NA
          method<-'HMR'
          advarsel<-'No residual degrees of freedom with n=3'
        }
      } else
      {
        if (n>3)
        {
          s2<-sum((konc-phi-f0.est*x)*(konc-phi-f0.est*x))/(n-3)
          y<-(kappa*tid+1)/(-kappa)
          O11<-n
          O12<-sum(x); O21<-O12
          O13<-sum(x*y); O31<-O13
          O22<-sum(x*x)
          O23<-sum(x*x*y); O32<-O23
          O33<-sum(x*x*y*y)
          O<-matrix(c(O11,O12,O13,O21,O22,O23,O31,O32,O33),nrow=3,ncol=3,byrow=TRUE)
          tjek<-try(solve(O),silent=TRUE)
          if (class(tjek)!='try-error')
          {
            f0.se<-sqrt(s2*solve(O)[2,2])
            f0seOK<-TRUE
          } else
          {
            f0.se<-NA
            f0seOK<-FALSE
          }
        } else {f0.se<-NA; f0seOK<-FALSE}
        if ((n>3)&(f0seOK))
        {
          f0.p<-2*pt(q=-abs(f0.est/f0.se),df=n-3)
          fraktil<-qt(p=0.975,df=n-3)
          f0.lo95<-f0.est-fraktil*f0.se
          f0.up95<-f0.est+fraktil*f0.se
          method<-'HMR'
          advarsel<-'None'
        } else
        {
          f0.p<-NA; f0.lo95<-NA; f0.up95<-NA
          method<-'HMR'
          if ((n<=3)&(!(f0seOK))) {advarsel<-'No residual degrees of freedom with n=3'}
          if ((n<=3)&(f0seOK))    {advarsel<-'No residual degrees of freedom with n=3'}
          if ((n>3)&(!(f0seOK)))  {advarsel<-'Standard error could not be computed'}
        }
      }
    } else
    {
      if (dacode==2)
      {
        # LR
        dum<-lsfit(tid/h,konc,intercept=TRUE)
        C0.est<-as.numeric(dum$coef[1])
        f0.est<-as.numeric(dum$coef[2])
        f0.se<-as.numeric(ls.diag(dum)$std.err)[2]
        f0.p<-2*pt(q=-abs(f0.est/f0.se),df=n-2)
        fraktil<-qt(p=0.975,df=n-2)
        f0.lo95<-f0.est-fraktil*f0.se
        f0.up95<-f0.est+fraktil*f0.se
        method<-'LR'
        if (C0.est<=0) {advarsel<-'C0<=0'} else {advarsel<-'None'}
      } else
      {
        # No flux
        f0.est<-0; f0.se<-NA; f0.p<-NA; f0.lo95<-NA; f0.up95<-NA; method<-'No flux'; advarsel<-'None'
      }
    }
  } else {f0.est<-NA; f0.se<-NA; f0.p<-NA; f0.lo95<-NA; f0.up95<-NA; method<-'None'; advarsel<-'Cancelled'}
  # Hvis "LR.always"
  if (LR.always)
  {
    if (dacode==2)
    {
      LR.f0<-f0.est
      LR.f0.se<-f0.se
      LR.f0.p<-f0.p
      LR.f0.lo95<-f0.lo95
      LR.f0.up95<-f0.up95
      LR.advarsel<-advarsel
    } else
    {
      dum<-lsfit(tid/h,konc,intercept=TRUE)
      LR.C0<-as.numeric(dum$coef[1])
      LR.f0<-as.numeric(dum$coef[2])
      LR.f0.se<-as.numeric(ls.diag(dum)$std.err)[2]
      LR.f0.p<-2*pt(q=-abs(LR.f0/LR.f0.se),df=n-2)
      LR.fraktil<-qt(p=0.975,df=n-2)
      LR.f0.lo95<-LR.f0-LR.fraktil*LR.f0.se
      LR.f0.up95<-LR.f0+LR.fraktil*LR.f0.se
      if (LR.C0<=0) {LR.advarsel<-'C0<=0'} else {LR.advarsel<-'None'}
    }
  } else {LR.f0<-NA; LR.f0.se<-NA; LR.f0.p<-NA; LR.f0.lo95<-NA; LR.f0.up95<-NA; LR.advarsel<-NA}
  
  ### Output
  list(serie=serie,f0=f0.est,f0.se=f0.se,f0.p=f0.p,f0.lo95=f0.lo95,f0.up95=f0.up95,PHMR.ptid=PHMR.ptid,PHMR.pkonc=PHMR.pkonc,
  method=method,warning=advarsel,LR.f0=LR.f0,LR.f0.se=LR.f0.se,LR.f0.p=LR.f0.p,LR.f0.lo95=LR.f0.lo95,LR.f0.up95=LR.f0.up95,LR.warning=LR.advarsel)
}
