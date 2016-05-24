
 GQD.plot=function(x,thin=1,burns,h=FALSE)
 {
    if(class(x)=='GQD.mcmc')
    {
      par.matrix=t(x$par.matrix)
      prop.matrix=t(x$prop.matrix)
      theta=par.matrix[,1]
      acc=x$acceptance.rate

      tt=1:length(acc)
      if(missing(burns)){burns =min(round(dim(x$par.matrix)[1]/2),25000)}

      tt = seq(burns,dim(x$par.matrix)[1],thin)

      #windows()
      nper=length(theta)
      if(nper==1){par(mfrow=c(1,2));d1=1;d2=1;}
      if(nper==2){par(mfrow=c(2,2));d1=1;d2=2;}
      if(nper==3){par(mfrow=c(2,2));d1=2;d2=2;}
      if(nper>3)
      {
        d1=1:((nper)+1)
        d2=d1
        O=outer(d1,d2)
        test=O-((nper)+1)
        test[test<0]=100
        test=test[1:4,1:4]
        test
        wh=which(test==min(test))
        
        d1=d1[col(test)[wh[1]]]
        d2=d2[row(test)[wh[1]]]
        par(mfrow=c(d1,d2))
      }
      cols=rainbow_hcl(nper, start = 10, end = 275,c=100,l=70)
      ylabs=paste0('theta[',1:nper,']')
      for(i in 1:nper)
      {
          if(h)
          {
            hist(par.matrix[i,tt],col=cols[i],main=ylabs[i],freq=F)
          }else
          {
           plot(prop.matrix[i,tt]~tt,col='gray90',type='s',main=ylabs[i],xlab='Iteration',ylab='')
           lines(par.matrix[i,tt]~tt,col=cols[i],type='s')
           abline(v=1,lty='dotdash')
          }
      }
      plot(acc,type='l',ylim=c(0,1),col='darkblue',main='Accept. Rate',xlab='Iteration',ylab='%/100')
      abline(h=seq(0,1,1/10),lty='dotted')
      #abline(v=burns,lty='dotdash')
      abline(h=0.4,lty='solid',col='red',lwd=1.2)
      abline(h=0.2,lty='solid',col='red',lwd=1.2)

      box()
   }
   
    if(class(x)=='GQD.mle')
   {
      stop('Nothing to plot!')
   }
   
   if(class(x)=='GQD.density')
   {
       if(requireNamespace('rgl', quietly = TRUE))
      {
        open3d(windowRect=c(50,50,640+50,50+640),zoom=0.95)
        persp3d(x=x$Xt,y=x$time,z=x$density,col=3,box=F,xlab='State (X_t)',ylab='Time (t)',zlab='Density f(X_t|X_s)')
        play3d(spin3d(axis=c(0,0,1), rpm=3), duration=10)
      }else
      {
        persp(x=x$Xt,y=x$time,z=x$density,col=3,xlab='State (X_t)',ylab='Time (t)',zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)
      }
   }
   

   if(class(x)=='BiGQD.density')
   {
     #x11()
     for(i in 1:dim(x$density)[3])
    {
      # Now illustrate the density:
      filled.contour(x$Xt,x$Yt,x$density[,,i],
      main=paste0('Transition Density \n (t = ',x$time[i],')'),
      color.palette=colorRampPalette(c('white','green','blue','red'))
      ,xlab='Xt',ylab='Yt')
    }
   }
 }