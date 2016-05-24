`PlotPTsmooth` <-
function(paz, pdip, x=0, y=0, siz=1, bcol='white', border="black", IMAGE=TRUE, CONT=TRUE, cont.col = 'black', pal=terrain.colors(100), LABS=FALSE, add=FALSE, NCP=50, NIP=200 )
  {
    
    ############  plot a smoothed density plot of the P and T axes
    ##########       smoothed focal mechanism summary
    if(missing(LABS)) { LABS=FALSE }
    
    if(missing(x)) { x=0 }
    if(missing(y)) { y=0 }
    if(missing(siz)) { siz=1 }
    if(missing(add)) { add=FALSE }
    if(missing(bcol)) { bcol='white' }
    if(missing(border)) { border='black' }
    if(missing(pal)) { pal=terrain.colors(100) }
    if(missing(IMAGE)) { IMAGE=TRUE }
    if(missing(CONT)) { CONT=FALSE }
    if(missing(cont.col)) {cont.col = 'black'}
       if(missing(NCP)) { NCP = 50 }
        if(missing(NIP)) { NIP=200 }

    
    PZZ =  focpoint(paz, pdip, col='red',  pch=3, lab="", UP=FALSE, PLOT=FALSE)
    KP = MASS::kde2d(PZZ$x, PZZ$y, n=NCP, lims=c(-1, 1, -1, 1))

    
    CC = PLTcirc(PLOT=FALSE, add=TRUE,  ndiv=72,  angs=c(-pi, pi))

    dx = KP$x[2] - KP$x[1]
    dy = KP$y[2] - KP$y[1]
    blankx = c(x+siz*CC$x, x+siz*CC$x[1], x+siz*(-1)-(dx), x+siz*(1)+(dx), x+siz*(1)+(dx), x+siz*(-1)-(dx), x+siz*(-1)-(dx))
    blanky = c(y+siz*CC$y, y+siz*CC$y[1], y+siz*(-1)-(dy), y+siz*(-1)-(dy), y+siz*(1)+(dy), y+siz*(1)+(dy), y+siz*(-1)-(dy))

    if(IMAGE==TRUE)
      {
       ##### polygon(x+siz*CC$x,  y+siz*CC$y , border=bcol,col=bcol)
       
         KP = MASS::kde2d(PZZ$x, PZZ$y, n=NIP, lims=c(-1, 1, -1, 1))
         M  = RPMG::meshgrid(KP$x, KP$y)
         flag = sqrt(M$x^2+M$y^2)>1
         KP$z[flag] = NA

    image(x+siz*KP$x,  y+siz*KP$y, KP$z, add=TRUE, col=pal, xpd=TRUE)
    
    ##### polygon(blankx, blanky , border=bcol,col=bcol)

  }

    if(CONT==TRUE)
      {
        if(IMAGE==FALSE & add==FALSE )  polygon(x+siz*CC$x,  y+siz*CC$y , border=bcol,col=bcol, xpd=TRUE)
              #####  polygon(x+siz*CC$x,  y+siz*CC$y , border=border,col=bcol)

         
        ##  contour(x+siz*KP$x,  y+siz*KP$y, KP$z, add=TRUE, col=cont.col, drawlabels =FALSE)

         cline.list=contourLines(x+siz*KP$x,  y+siz*KP$y, KP$z)
         templines <- function(clines) {
           cx = clines[[2]]
           cy = clines[[3]]
           flag = sqrt((cx-x)^2+(cy-y)^2)> siz
           cx[flag] = NA
           cy[flag] = NA
           lines(cx, cy, col=cont.col, xpd=TRUE)
         }
         
         invisible(lapply(cline.list, templines))
         
         
         
       }

    if(!is.na(border)) { lines(x+siz*CC$x,  y+siz*CC$y , col=border, xpd=TRUE) }
    if(FALSE)
      {

        focpoint(paz, pdip, col='red',  pch=3, lab="", UP=FALSE)
###  title text(0,1.04,labels="P-axes 2D Density", font=2, cex=1.2)

      }
    
    if(LABS==TRUE)
      {
        
      }
    
    
  }

