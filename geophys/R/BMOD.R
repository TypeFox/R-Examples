BMOD<-function(bill, nstn=100, PLOT=TRUE, obs=NULL )
  {
#######  two parts - model and plotting 
    if(missing(nstn)) nstn = 100
    if(missing(PLOT)) PLOT=TRUE

    ZCOLS = RPMG::pastel.colors(12, seed=1 )

    DZ = abs( bill$zmax - bill$zmin)
    
    PZ = DZ*0.2
    
    ZMAX = bill$zmin+PZ
    ZMIN = bill$zmin+DZ*0.02
    
    xstart = bill$xmin
    xend =   bill$xmax

    xs = seq(from=xstart, by=(xend-xstart)/nstn , length=nstn)
    zs = rep(0, length=length(xs))
    
    den = 0.2
    n = length(bill$mod)
    
    ALLgravZ = rep(0, length(xs))

    for(i in 1:bill$n)
      {
        den = bill$mod[[i]]$rho
        if(is.na(den)) den=0.2
        if(is.null(den)) den=0.2

        
        pol = flipZEE(bill$mod[[i]])

        dc = dircheck(pol)
        ## print(dc)
        L1 = length(which(dc>0))
        if(L1>=(length(pol$x)/2))
          {
            cat(paste("Reversing", i, sep=":") , sep="\n")
            pol =  rev2RH(pol)
          }
        
        Ngrav = DGzx(xs, zs,  pol$x, pol$y, den)
        ALLgravZ = ALLgravZ + Ngrav$Gz
      }
    ##  plot(xs, ALLgravZ )
#########################
    if(PLOT==TRUE)
      {
        plot(c(bill$xmin, bill$xmax), c( bill$zmax, ZMAX), type='n' , axes=FALSE, ann=FALSE )
        grid()
        title(xlab="X-m", ylab="Depth, m")
        rect(bill$xmin, bill$zmin, bill$xmax, bill$zmax,  col=NA, border="black", xpd=TRUE )

        for(i in 1:bill$n)
          {

            den = bill$mod[[i]]$rho
            pol = (bill$mod[[i]])
            polygon(pol$x, pol$y, col=ZCOLS[i], border="black")

###  need to check to make sure the polygon is right handed
            
            text(pol$x, pol$y, labels=1:length(pol$x))
            cenP =  centroid(pol)
            text( cenP[1], cenP[2], paste(i,":", den, sep="" ))
          }

        g1 = min(ALLgravZ)
        g2 = max(ALLgravZ)
        RGz = RESCALE(ALLgravZ, ZMIN, ZMAX, g1, g2)

        rect( bill$xmin, ZMIN, bill$xmax, ZMAX, col='white', border=grey(.9) )

        axis(1)
        paz = pretty(c( bill$zmin, bill$zmax))
        axis(2, at=paz)

        ## gravaxis
        gzmin = min(ALLgravZ)
        gzmax = max(ALLgravZ) 

        text(bill$xmax, ZMIN, labels=format(gzmin, digits=4), pos=4, xpd=TRUE )
        text(bill$xmax, ZMAX, labels=format(gzmax, digits=4), pos=4, xpd=TRUE )

        pgz = pretty(range(ALLgravZ))
        pgz = pgz[ pgz>gzmin &  pgz<gzmax  ] 
        rpgz = RESCALE(pgz, ZMIN, ZMAX, g1, g2)

        segments(bill$xmin, rpgz, bill$xmax, rpgz, lty=2, col=grey(.9) )
        lines(xs, RGz, col='blue')


        if(!is.null(obs))
          {
               ## obsGz = RESCALE(obs$g, ZMIN, ZMAX, g1, g2)
               obsGz = RESCALE(obs$g, ZMIN, ZMAX, min(obs$g), max(obs$g) )
               points(obs$x, obsGz, pch=3, col='red')
               text(bill$xmin, ZMIN, labels=format(min(obs$g), digits=4), pos=2, xpd=TRUE, col='red' )
               text(bill$xmin, ZMAX, labels=format( max(obs$g), digits=4), pos=2, xpd=TRUE, col='red' )

          }

        
      }

    invisible(ALLgravZ)
  }
