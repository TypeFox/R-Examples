### Function for plotting GLCM-matrices

plotGLCM<-function(GLCM)
{
plot(0,0, xlim=c(0,255), ylim=c(0,255), type="n", main=paste("GLCM-Matrix"))
for(x in 1:255)
{
   for(y in 1:255)
   {
      if(GLCM[x,y]>0)
      { 
          points(x,y, cex=0.25, pch=19)
      }    
    }  
}

}##eo function plotGLCM
