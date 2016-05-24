figureEF <-
function(boundary, K, sm=3, lg=6)
{
  par(ask=T)
  for(i in 1:K)
  {             
    plot(0,0, xlim=c(-sm,lg), ylim=c(-lg,sm), xlab="",ylab="",col='white',bty="n",axes=F)
    mtext(paste("stage ",as.character(i)),1,0)
    
    polygon(c(-lg, lg,lg), c(-lg, -lg, lg), border = NA,
            col=rgb(0.95,0.8,0.8));   
    
    polygon(c(boundary$futilL[i],lg,lg, boundary$futilL[i]),
            c(boundary$futilU[i],boundary$futilU[i], -lg,-lg),
            border = NA,col=rgb(0.9,0.9,0.9))   
    
    polygon(c(boundary$equivL[i],lg,lg, boundary$equivL[i]),
            c(boundary$equivU[i],boundary$equivU[i], -lg,-lg),
            border = NA,col=rgb(0.8,0.9,0.8))     
    
    segments(boundary$equivL[i],boundary$equivU[i], lg, boundary$equivU[i], 
             lty=2, lwd=2, col='green4')
    segments(boundary$equivL[i],boundary$equivU[i], boundary$equivL[i], -lg, 
             lty=2, lwd=2, col='green4')    
    points(boundary$equivL[1:i],boundary$equivU[1:i], pch=8)
    
    if(i<K){
      segments(boundary$futilL[i],boundary$futilU[i], boundary$futilL[i], -lg, 
               lty=2, lwd=2, col='red')
      segments(boundary$futilL[i],boundary$futilU[i], lg, boundary$futilU[i], 
               col='red', lty=2, lwd=2)     
    }
    polygon(c(-sm,-sm,sm),c(-sm,sm,sm), col='white', border = NA)   
    
    points(boundary$equivL[1:i],boundary$equivU[1:i], pch=8, col='green4')
    points(boundary$futilL[1:i],boundary$futilU[1:i], pch=1, col='red')
    legend(-sm,sm, pch=c(-1,8,1),col=c('white','green4','red'),
           bg="transparent", bty = "n",
            c('critical values for','equivalence','futility'))
    
    text(lg*2/3,-lg*2/3,"equivalence") 
    text(lg*2/3,sm/2,"futility")   
    if(i<K) text(lg*2/3,boundary$equivU[i]/2,"continuation")
   
    segments(-sm,0,lg,0); text(lg,0.5,"t(L)")
    segments(0,-lg,0,sm); text(-0.5,sm,"t(U)")
    abline(0,1, lty=3,col='white')
  }        
}
