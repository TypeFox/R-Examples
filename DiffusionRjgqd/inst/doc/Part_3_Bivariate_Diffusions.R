## ----fig.align='center'--------------------------------------------------

 library(DiffusionRjgqd)
 JGQD.remove()
 a00 = function(t){0.5*5}
 a10 = function(t){-0.5}
 a01 = function(t){-0.1}
 c10 = function(t){0.4^2}

 b00 = function(t){0.4*6}
 b01 = function(t){-0.4}
 b10 = function(t){-0.1}
 f01 = function(t){0.3^2}

 Lam00= function(t){1}

 Jmu1=function(t){0.75*(1+sin(2*pi*t))}
 Jmu2=function(t){0.75*(1+sin(2*pi*t))}
 Jsig11=function(t){0.75^2*(1+0.8*sin(2*pi*t))^2}
 Jsig22=function(t){0.75^2*(1+0.8*sin(2*pi*t))^2}

 xx=seq(3,11,1/10)
 yy=seq(3,11,1/10)
 res= BiJGQD.density(7,7,xx,yy,0,1,1/100,Dtype='Saddlepoint')

## ----fig.align = 'center'------------------------------------------------
 #Simulate the process
 mux     = function(x,y,t){a00(t)+a10(t)*x+a01(t)*y}
 sigmax  = function(x,y,t){sqrt(c10(t)*x)}
 muy     = function(x,y,t){b00(t)+b10(t)*x+b01(t)*y}
 sigmay  = function(x,y,t){sqrt(f01(t)*y)}
 lambda1 = function(x,y,t){Lam00(t)}
 lambda2 = function(x,y,t){rep(0,length(x))}

 j11     = function(x,y,z){z}
 j12     = function(x,y,z){z}
 j21     = function(x,y,z){z}
 j22     = function(x,y,z){z}
 simulate=function(x0=7,y0=7,N=10000,TT=5,delta=1/1000,pts,brks=30,plt=FALSE)
 {
   library(colorspace)
   colpal=function(n){rev(sequential_hcl(n,power=0.8,l=c(40,100)))}

   d=0                  # Time index
   tt=seq(0,TT,delta)   # Time sequance
   X=rep(x0,N)          # Initialize state vectors
   Y=rep(y0,N)



   x.traj = rep(x0,length(tt))
   y.traj = rep(y0,length(tt))
   x.jump = rep(0,length(tt))
   y.jump = rep(0,length(tt))

   # Storage for histogram snapshots:
   evts = rep(0,N)
   for(i in 2:length(tt))
   {

    X=X+mux(X,Y,d)*delta+sigmax(X,Y,d)*rnorm(N,sd=sqrt(delta))
    Y=Y+muy(X,Y,d)*delta+sigmay(X,Y,d)*rnorm(N,sd=sqrt(delta))
    events1 = (lambda1(X,Y,d)*delta>runif(N))
    if(any(events1))
    {
      wh=which(events1)
      evts[wh]=evts[wh]+1
      X[wh]=X[wh]+j11(X[wh],Y[wh],rnorm(length(wh),Jmu1(d),sqrt(Jsig11(d))))
      Y[wh]=Y[wh]+j21(X[wh],Y[wh],rnorm(length(wh),Jmu2(d),sqrt(Jsig22(d))))
    }
    events2 = (lambda2(X,Y,d)*delta>runif(N))

    d=d+delta
    if(sum(round(pts,3)==round(d,3))!=0)
    {
    if(plt)
    {
     expr1 = expression(X_t)
     expr2 = expression(Y_t)
     color.palette=colorRampPalette(c('green','blue','red'))
     filled.contour(res$Xt,res$Yt,res$density[,,i],
                    main=paste0('Transition Density \n (t = ',round(d,2),')'),
                    color.palette=colpal,
                    nlevels=41,xlab=expression(X[t]),ylab=expression(Y[t]),plot.axes=
     {
        # Add simulated trajectories
        points(Y~X,pch=c(20,3)[(evts>0)+1],col=c('black','red')[(evts>0)+1],cex=c(0.9,0.6)[(evts>0)+1])
        if(any(events2))
        {
          wh=which(events2)
          segments(xpreee[wh],ypreee[wh],X[wh],Y[wh],col='gray')
        }
        axis(1);axis(2);
        # Add a legend
        legend('topright',col=c('black','red'),pch=c(20,3),
                legend=c('Simulated Trajectories','Jumped'))
        yy=contourLines(res$Xt,res$Yt,res$density[,,i],levels=seq(0.01,0.1,length=10))
        if(length(yy)>0)
        {
        for(j in 1:length(yy))
        {
         lines(yy[[j]])
        }
        }
     })
    }
    }
  }
}
sim=simulate(7,7,N=200,TT=0.75,delta=1/100,plt=TRUE,pts=c(0.13,0.28,0.38,0.51,0.63,0.75))


## ----eval=FALSE----------------------------------------------------------
#  browseVignettes('DiffusionRjgqd')

