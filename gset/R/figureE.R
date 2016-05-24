figureE <-
function(boundary, K, sm=3, lg=6)
{
  par(ask=T)
  for(i in 1:K)
  {             
    plot(0,0, xlim=c(-sm,lg), ylim=c(-lg,sm), xlab="",ylab="",
         col='white',bty="n",axes=F)
    mtext(paste("stage ",as.character(i)),1,0)
    
    polygon(c(boundary$equivL[i],lg,lg, boundary$equivL[i]),
            c(boundary$equivU[i],boundary$equivU[i], -lg,-lg),
            border = NA,col=rgb(0.8,0.9,0.8))     

    polygon(c(boundary$equivL[i],-sm,-sm, lg, lg, boundary$equivL[i]),
            c(-lg,-lg,sm,sm,boundary$equivU[i],boundary$equivU[i]),
            border = NA,col=rgb(0.9,0.9,0.9)) 
    
    polygon(c(-sm,-sm,sm),c(-sm,sm,sm), col='white', border = NA)   
    
    segments(boundary$equivL[i],boundary$equivU[i], lg, boundary$equivU[i], 
             lty=2, lwd=2, col='green4')     
    segments(boundary$equivL[i],boundary$equivU[i], boundary$equivL[i],-lg, 
             lty=2, lwd=2, col='green4')       
    
    segments(-sm,0,lg,0); text(lg,0.5,"t(L)")
    segments(0,-lg,0,sm); text(-0.5,sm,"t(U)")
    
    points(boundary$equivL[1:i],boundary$equivU[1:i], pch=8, col='green4')

    legend(-sm,sm, pch=c(8, -1),col=c('green4','white'),
           bg="transparent", bty = "n",
           c('critical values for','equivalence'))
    text(lg*2/3,-lg*2/3,"equivalence") 
    if(i<K)  text(lg*2/3,sm/2,"continuation")
    else text(lg*2/3,sm/2,"non-equivalence")
    
    abline(0,1, lty=3,col='white')
  }        
}
