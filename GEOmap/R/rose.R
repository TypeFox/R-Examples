
rose<- function(angles, bins, x=0, y=0, col="black", border="black", annot=FALSE, main = "", prop = 1, pts = FALSE,
                cex = 1, pch = 16, dotsep = 40, siz = 1, LABS=LABS, LABangle=180, add=FALSE, SYM=FALSE)
{
###  angles = angles
###  bins = 

  if(missing(bins)) bins =36
    if(missing(x)) x=0
    if(missing(y))  y=0
  if(missing(col)) col="black"
  if(missing(border)) border="black"
  if(missing(main)) main = "" 
  if(missing(prop)) prop = 1
  if(missing(pts))  pts = FALSE
  if(missing(cex))  cex = 1
  if(missing(pch)) pch = 16
  if(missing(dotsep)) dotsep = 40
  if(missing(siz)) siz = 1
  if(missing(annot)) annot=TRUE
  if(missing(LABangle)) LABangle=NULL
  if(missing(add)) add=FALSE
   if(missing(SYM)) SYM=FALSE
  

  if(missing(LABS)) { LABS = c("90", "270", "180", "0") }
                                        #   LABS = c("N", "S", "W", "E")

fmod<-function(k, m)
  {
    j = floor(k/m)
    a = k-m*j
    return(a)
  }

  ###if(SYM)
  ###  {
   ###   angles =   c(angles, angles+pi)

 ###   }
  
  angles <- RPMG::fmod(angles , (2 * pi))

  cx <- cos(seq(0, 2 * pi, length = 1000))
  sy <- sin(seq(0, 2 * pi, length = 1000))

  if(add==FALSE)
    {
     ###  plot(cx,sy,
       ###     axes = FALSE, xlab = "", ylab = "",
      ###      main = main, type = "n", asp=1, xlim = siz * c(-1, 1), ylim = siz* c(-1, 1))
      par(mar=c(0,0,0,0))
      
      plot(siz*cx, siz*sy,
           axes = FALSE, xlab = "", ylab = "",
           main = main, type = "n", asp=1)

    }

  up <- par("usr")
  ui <- par("pin")
  
  ratx  <-(up[2]-up[1])/ui[1]
  raty <- (up[4]-up[3])/ui[2]
  
  usizx  <-2*siz*ratx
  usizy  <-2*siz*raty
  
  if(annot==TRUE)
    {
      lines(x+usizx*cx,y+usizy*sy)
      
      labplacement <- 0.95
      
      lines(x+usizx*c(0, 0), y+usizy*c(labplacement, 1))
      text(x+usizx*0.0,y+usizy* labplacement, LABS[1], pos=1, cex = 1.1)
      lines(x+usizx*c(0, 0), y+usizy*c(-labplacement, -1))
      text(x+usizx*0.0, y+usizy*(-labplacement), LABS[2], pos=3, cex = 1.1)
      lines(x+usizx*c(-1, -labplacement), y+usizy*c(0, 0))
      text(x+usizx*(-labplacement), y+usizy*0, LABS[3], pos=4, cex = 1.1)
      lines(x+usizx*c(labplacement, 1), y+usizy*c(0, 0))
      text(x+usizx*labplacement, y+usizy*0, LABS[4], pos=2, cex = 1.1)

      
    }
  ##  points(0, 0, cex = 1)
  n <- length(angles)
  freq <- rep(0, bins)
  arc <- (2 * pi)/bins
  darcs = (0:bins)*arc


  if(SYM)
    {

      jx = cos(angles)
      jy = sin(angles)
      
      angles = atan2(abs(jy), jx)
      

      

      ###  flip the vectors
      ihalf = findInterval(pi, darcs)
      
      i = 1
      freq[i] <- sum(  angles <= darcs[i+1] & angles >= darcs[i]  )

      jind = 2:ihalf
      for(i in  jind ) {
        freq[i] <- sum( angles <= darcs[i+1] & angles > darcs[i]   )
      }

      
      freq[(ihalf+1):bins] = (freq[1:ihalf])
      

    }
  else
    {

      i = 1
        freq[i] <- sum(  angles <= darcs[i+1] & angles >= darcs[i]  )
      for(i in 2:bins) {
        freq[i] <- sum(  angles <= darcs[i+1] & angles > darcs[i]  )
      }


    }






  
  rel.freq <- freq/n

  wmax <-	which.max(rel.freq)

  prop <- 1/sqrt(rel.freq[wmax])	

  maxnum <- freq[wmax]

  pret <-	pretty(1:maxnum)
  pret <- pret[pret>0 & pret<=maxnum]

  if(!is.null(LABangle)) LABangle  <- pi*LABangle/180
  if(annot==TRUE)
    {
      for(i in 1:length(pret))
        {
          r <- sqrt(pret[i]/n) * prop
          lines(x+usizx*r*cx, y+usizy*r*sy, col=grey(.8))
          if(!is.null(LABangle))  text(x+usizx*r*cos(LABangle), y+usizy*r*sin(LABangle), labels=pret[i])
        }
      
    }
#########   for area scaling (correct) use this:
  
  radius <- sqrt(rel.freq) * prop
  
##########   for non feq proportional to length use this:
  ##  radius <- (rel.freq) * prop

  sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
  mids <- seq(arc/2, 2 * pi - pi/bins, length = bins)
  index <- cex/dotsep
  for(i in 1:bins) {
    if(rel.freq[i] != 0) {
      
      pieslice<- list(x= x+usizx*c(0, radius[i] * cos(sector[i]), radius[i] * cos(sector[i]), radius[i] * cos(sector[i] + (2 * pi)/bins), 0)  ,
                      y= y+usizy*c(0, radius[i] * sin(sector[i]) , radius[i] * sin(sector[i]), radius[i] * sin(sector[i] + (2 * pi)/bins),  0 )   )


      polygon(pieslice, col=col, border=border)

      
      if(pts == TRUE) {
        for(j in 0:(freq[i] - 1)) {
          r <- 1 + j * index
          px <- x+usizx*r * cos(mids[i])
          py <- y+usizy*r * sin(mids[i])
          points(px, py, cex = cex, pch = pch)
        }
      }
    }
  }
  return(list(sector=sector, radius=radius, usizx=usizx, usizy=usizy, x=x, y=y ))
}

