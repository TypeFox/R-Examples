plotJACKLLZ<-function(hjack, sta, proj=NULL, PLOT=0, PS=FALSE, fbase="jack", width = c(10, 5), height = c(8,8)  )
{
###  there are 2 possible plot created eath with a different default size

  ###  if plot = {0,1,2}
#######  plot the output of HiJACK
  JPNG<-function(file="tmp", width = 8, height = 8)
    {
      
      png(filename = file,
          width = width, height = height, units='in', pointsize = 12, res=300)
    }

  if(length(width)<2) { width = c(width[1], 0.5*width[1]) }
  if(length(height)<2) { height = c(height[1], height[1]) }

  
  YEYE=hjack$Y
  XEYE=hjack$X
  ZEYE=hjack$Z
  
  L1 = as.vector( unlist( lapply(YEYE, "length") ) )
  Lw = which(L1>0)

  YEYEB = vector(mode="list")
  for(i in 1:length(Lw) )
    {
      YEYEB[[i]] = YEYE[[Lw[i]]]

    }
  names( YEYEB) = names(YEYE)[Lw]


 
################################

  L1 = as.vector( unlist( lapply(XEYE, "length") ) )
  Lw = which(L1>0)

  XEYEB = vector(mode="list")
  for(i in 1:length(Lw) )
    {
      XEYEB[[i]] = XEYE[[Lw[i]]]

    }
  names( XEYEB) = names(XEYE)[Lw]
################################
  L1 = as.vector( unlist( lapply(ZEYE, "length") ) )
  Lw = which(L1>0)

  ZEYEB = vector(mode="list")
  for(i in 1:length(Lw) )
    {
      ZEYEB[[i]] = ZEYE[[Lw[i]]]

    }
  names( ZEYEB) = names(ZEYE)[Lw]


##  boxplot.stats(x, coef = 1.5, do.conf = TRUE, do.out = TRUE)
Blat = boxplot(YEYEB, plot=FALSE)
  
  Blon = boxplot(XEYEB, plot=FALSE)
 
  Bz = boxplot(ZEYEB, plot=FALSE)

  if(PS==FALSE)
    {
  fn1=NA
  fn2=NA
}
 ##
  if( PLOT==0 | PLOT==1)
    {
      
      if(PS==TRUE)
        {
          fn1=paste(fbase,"1.png", sep="" )
          JPNG(file=fn1, width=width[1], height=height[1] )
        }
      else
        {
          op <- par(no.readonly = TRUE) 
          dev.new(width=width[1], height=height[1])
        }

      
      par(mfrow=c(3,1))

      
      bxp(Blat, varwidth=TRUE)
      title("Station Influence Lat")

      bxp(Blon, varwidth=TRUE)
      title("Station Influence Lon")

      bxp(Bz, varwidth=TRUE)
      title("Station Influence Z")
      if(PS==TRUE) { dev.off() }
      else
        {
          par(op)
        }
    }


  if( PLOT==0 | PLOT==2)
    {

      if(PS==TRUE)
        {
          fn2=paste(fbase,"2.png", sep="" )
          JPNG(file=fn2, width=width[2], height=height[2] )
        }
      else
        {
          dev.new(width=width[2], height=height[2])
          op <- par(no.readonly = TRUE) 
        }




  par(mfrow=c(3,1))
     

  Gxy = GEOmap::GLOB.XY(  sta$lat   ,sta$lon , proj)
  plot(Gxy$x, Gxy$y, type='n', xlab="km", ylab="km", asp=1 )
  imageINFLUENCE(Blat, sta, proj)
  title("Station Influence Latitude")
 
  plot(Gxy$x, Gxy$y, type='n', xlab="km", ylab="km", asp=1 )
  imageINFLUENCE(Blon, sta, proj)
  title("Station Influence Longitude")
 
  plot(Gxy$x, Gxy$y, type='n', xlab="km", ylab="km", asp=1 )
  imageINFLUENCE(Bz, sta, proj)
  title("Station Influence Depth")


      if(PS==TRUE) { dev.off()}
      else
        {
          par(op)
        }
      
    }
  return(list( X=XEYEB, Y=YEYEB, Z=ZEYEB, files=c(fn1, fn2)  ) )
}
