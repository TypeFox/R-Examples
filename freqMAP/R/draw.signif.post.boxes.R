`draw.signif.post.boxes` <-
function(post.dat,ylims,alpha1,alpha2,color1,color2,
                                   dens1,dens2,ylog=FALSE){

  if(alpha2>alpha1|alpha1>.5|alpha2>.5)
    stop("invalid: alpha1 must be g.t. alpha2")

  #sort post.dat by the first column
  post.dat <- post.dat[order(post.dat[,1]),]
  
  #plot the alpha2 (stronger p-value) boxes second
  for(alpha in c(alpha1,alpha2)){
    for(i in 1:nrow(post.dat)){
      if(!is.na(post.dat[i,2]) & (post.dat[i,2]>(1-alpha/2)|post.dat[i,2]<(alpha/2))){
        
        boxheight <- (ylims[2]-ylims[1])/20
        if(ylog) boxheight=exp((log(ylims[2])-log(ylims[1]))/20)
        
        if(alpha==alpha1){
          signif.color <- color1
          signif.dens <- dens1
        }
        if(alpha==alpha2){
          signif.color <- color2
          signif.dens <- dens2
        }
        y<-c(ylims[1],ylims[1],ylims[1]+boxheight,ylims[1]+boxheight)
        if(ylog) y<-c(ylims[1],ylims[1],ylims[1]*boxheight,ylims[1]*boxheight)
        if(i==1){
          #far left of plot
          boxdims <- c(0,
                       (post.dat[i+1,1]-post.dat[i,1])/2,
                       (post.dat[i+1,1]-post.dat[i,1])/2,
                       0)
        }else if(i==nrow(post.dat)){
          #far right of plot
          boxdims <- c((post.dat[i-1,1]-post.dat[i,1])/2,
                       0,
                       0,
                       (post.dat[i-1,1]-post.dat[i,1])/2)
        }else{
          #somewhere in the middle
          boxdims <- c((post.dat[i-1,1]-post.dat[i,1])/2,
                       (post.dat[i+1,1]-post.dat[i,1])/2,
                       (post.dat[i+1,1]-post.dat[i,1])/2,
                       (post.dat[i-1,1]-post.dat[i,1])/2)
        }
        x<-post.dat[i,1]+boxdims
        #if on the left or right edge, then clip the box inside the plot range
        x[1] <- max(x[1],min(post.dat[,1],na.rm=TRUE))
        x[2] <- min(x[2],max(post.dat[,1],na.rm=TRUE))
        x[3] <- x[2]
        x[4] <- x[1]
        polygon(x=x,y=y,
                density=signif.dens,col=signif.color,angle=90,border=NA)
      }
    }
  }

}

