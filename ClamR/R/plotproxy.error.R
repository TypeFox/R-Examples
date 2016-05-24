plotproxy.error<-function(x,y,gout, type=1, xlim=NULL, ylim=NULL, ylab="", xlab="", main="" )
  {

###########  use output of proxyJK to plotall sinusoids
    jout = gout$JOUT
    omids=gout$omids
    pmids=gout$pmids
    delw  = gout$delw

    if(missing(ylab)) ylab = expression(delta*"18O(% VPDB)")
    if(missing(xlab)) xlab = "Distance from Margin (mm)"
    if(missing(main)) main=""

    if(missing(xlim) | is.null(xlim) )   xlim = range(x)
       if(missing(ylim) | is.null(ylim)  ) ylim = range(y)
  

    if(type==1)
      {
        plot(x,y, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, pch=15, col='blue')
        points(omids,  pmids, col='brown', pch=20)
        lines(omids,  pmids, col='brown', lwd=2)

        error.bar(omids,  pmids,  pmids-delw,   pmids+delw , pch = 20, col = 'brown', barw = 0.05 , add=TRUE)

      }
    if(type==2)
      {

        plot(range(c(x, omids) )  , range(c(y, pmids-delw, pmids+delw)) ,type='n',  xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, pch=15, col='blue')


        polygon(c(omids, rev(omids)) , c(pmids-delw, rev(pmids+delw)), col=grey(.8), border=NA)

        points(x,y, pch=15, col='blue')

        points(omids,  pmids, col='brown', pch=20)
        lines(omids,  pmids, col='brown', lwd=2)





      }
    

  }


