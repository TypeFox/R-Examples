slideshow<-function(P = c("hi", "there", "sugar pie" ), dy=0.2, EX=0.1, ht=3, font=2, anim=FALSE   )
  {

    
    N = length(P)

    if(missing(dy) ) dy = 1/N
    if(missing(EX) )   EX = 1/10
    if(missing(ht) )    ht = 3
    if(missing(font) )  font=2
    if(missing(anim) )  anim=FALSE 

    
    plot(c(0,1),c(0,1), type='n', ann=FALSE, axes=FALSE)

    op <- options();
    options(locatorBell=FALSE)
    
    for(i in 1:N)
      {
        why = 1-( (i-1) *dy)
        text(x=EX, y=why, labels=P[i], cex=ht, font=font, xpd=TRUE, pos=4)
        if(anim) locator(1)
      }
    options(op)     # reset (all) initial options

    
  }
