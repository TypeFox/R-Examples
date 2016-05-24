`plotzoeppritz` <-
function(A, zoepcols = c("red", "green" , "blue", "purple"), zoeplty=c(1,1,1,1)  )
  {

    if(missing(zoepcols) ) zoepcols = c("red", "green" , "blue", "purple")
    if(missing(zoeplty) ) zoeplty=c(1,1,1,1)

    #########  A = list(angle=angle, rmat=rmat, rra=rra, ria=ria, ang=ang)
    plot(range(A$angle), range(A$rmat, na.rm=TRUE), col=4, type="n", xlab="Angle of Incidence",
         ylab=A$chtype ) ;

    u = par("usr")
    
    A$rmat[ A$rmat==0.0 ] = NA
    
    for(i in 1:length(A$rmat[1,]))
      {
        
        lines(A$angle, A$rmat[,i], col=zoepcols[i], lty=zoeplty[i])
      }
    
    if(!is.null(A$alphacrit))  { abline(v=A$alphacrit, lty=2, col=grey(.75) ) }
    ##  title(paste(sep="", chincw, " Incident/", choutkind, " Out" ) );
    
    ## text(0 , 1.15, paste('Layer 1: Vp=',alpha1,' Vs=',beta1, ' Rho=', rho1), pos=4)
    ## text(0, 1.05, paste('Layer 2: Vp=',alpha2,' Vs=',beta2, ' Rho=', rho2), pos=4)
    
    ##  legend("topleft",   legend=c(paste(sep=" ", "Incident:", A$incw),
   ##                       paste('Layer 1:', bquote(alpha[1]), '=',alpha1, bquote(beta[1]), '=',beta1, bquote(rho[1]), '=', rho1),
    ##                      paste('Layer 2: Vp=',alpha2,' Vs=',beta2, ' Rho=', rho2)))

    
    
    leg1 = bquote(alpha[1] == .(alpha1))
    
    vskip = strheight(leg1, units = "user", cex = 1)*1.4
    
    text(u[1], u[4] -vskip, leg1, pos=4)
    leg1 = bquote(beta[1] == .(beta1))
    text(u[1], u[4] -2*vskip, leg1, pos=4)
    leg1 = bquote(rho[1] == .(rho1))
    text(u[1], u[4] -3*vskip, leg1, pos=4)
    
    leg1 = bquote(alpha[2] == .(alpha2))
    text(u[1], u[4] -4*vskip, leg1, pos=4)
    leg1 = bquote(beta[2] == .(beta2))
    text(u[1], u[4] -5*vskip, leg1, pos=4)
    leg1 = bquote(rho[2] == .(rho2))
    text(u[1], u[4] -6*vskip, leg1, pos=4)
    
    
    
    legend("topright",  legend=c("P-reflected", "S-reflected", "P-refracted", "S-refracted"), lty=zoeplty, lwd=2, col=zoepcols)

  }

