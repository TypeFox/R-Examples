varsquiggle<-function(GH, sel=c(1,2), WIN=c(0,1) , dist=NULL, thick=1 , FLIP=FALSE, filcol='blue', tracecol='blue')
  {
#####  plot a seismic section as in exploration seismology
####  with varian - squiggle display

    ##  source("varsquig.R")

    if(missing(sel)) { sel= 1:length(GH$JSTR) }

     dt1 = GH$dt[sel[1] ]

    if(missing(WIN)) {

      WIN = c(0,  length(GH$JSTR[[  sel[1] ]]   )*dt1    )


      }

    if(missing(filcol))
      {
        filcol= rep('blue', length(sel))


      }
    if(length(filcol)==1)
      {

        filcol= rep(filcol, length(sel))

      }

    if(missing(tracecol))
      {
        tracecol= rep('blue', length(sel))


      }
    if(length(tracecol)==1)
      {

        tracecol= rep(tracecol, length(sel))

      }

   if(missing(dist))
     {
       dist = seq(from=1, to=length(sel))


     }


    
 

    
    if(missing(thick))
      {
        thick = NULL
      }
   

    N = length(sel)

    XMAT  = NULL

    x = seq(from=0, length=length(GH$JSTR[[  sel[1] ]]   )   , by=dt1)

    
    j = 1
    for(i in 1:N)
      {
        sig = GH$JSTR[[  sel[i] ]]
        len = length(sig)
        x = seq(from=0, length=len , by=dt1)
      

        
        sig = sig[x>=WIN[1] & x<=WIN[2]]

        sig = sig-mean(sig)
        XMAT = cbind(XMAT, sig  )

      }


    d = dim(XMAT)
    x = seq(from=0, length=d[1], by=dt1)

 matsquiggle(XMAT, dt1, dist, thick=thick, FLIP=FALSE, filcol=filcol,tracecol=tracecol,  PLOT=FALSE, add=FALSE)
  grid(col=rgb(.8,.6,.6) , nx=20)
  axis(1)
  axis(2)
  u=par("usr")
  
  text(rep(u[2], times=length(dist)) , dist, labels=GH$STNS[sel], pos=2, xpd=TRUE)

  
  Msquig = matsquiggle(XMAT, dt1, dist, thick=thick, FLIP=FALSE, filcol=filcol,tracecol=tracecol,  add=TRUE, PLOT=TRUE)
  


    invisible( Msquig  )
  }
