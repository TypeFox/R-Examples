PostVquake<-function(MANYeq, GX, GY, XY, proj, add=FALSE, ... )
  {

    zcols = c('red', 'blue', 'purple', 'cyan' )

  
    ZEDS = matrix(ncol=3, nrow=length(MANYeq))
    

    for(i in 1:length(MANYeq))
      {
        zip = MANYeq[[i]]$Ksolutions
        ZEDS[i, ] = c(MANYeq[[i]]$EQ$lat,MANYeq[[i]]$EQ$lon, MANYeq[[i]]$EQ$z ) 
        
        for(iz in 1:length(zip))
          {
            miz =  zip[[iz]]
            points(miz[,1], miz[,2], col=zcols[iz], pch=iz) 
            lines(miz[,1], miz[,2], col=zcols[iz] )

          }

      }


    if(!add)
      {
    plot(GX, GY, type='n', xlab="km" , ylab="km" , asp=1)
    points(XY, pch=6)
  }
    ZEXY = GEOmap::GLOB.XY(ZEDS[,1], ZEDS[,2], proj)
    points(ZEXY, pch=8, col='red')



    for(i in 1:length(MANYeq))
      {
        KOV =  MANYeq[[i]]$ERR$cov[2:4, 2:4]
        ndf = MANYeq[[i]]$ERR$ndf
        
        eqlipse(ZEXY$x[i], ZEXY$y[i] , KOV,   wcols = c(1,2) , dof=ndf,  ...  )
      }
    
    

  }


######################################















