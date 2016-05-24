PLOTbox <-
function(Rax, Rbox, axcol= 'black', boxcol= 'blue' )
  {
    if(missing(axcol)){axcol= 'black'}
    if(missing(boxcol)){boxcol= 'blue'}

    
    segments(  Rax[1,1],Rax[1,2], Rax[2,1],Rax[2,2], col=axcol, xpd=TRUE)
    segments(  Rax[3,1],Rax[3,2], Rax[4,1],Rax[4,2], col=axcol, xpd=TRUE)
    segments(  Rax[5,1],Rax[5,2], Rax[6,1],Rax[6,2], col=axcol, xpd=TRUE)

###  labels for axes

    text(Rax[2,1],Rax[2,2], labels="x", pos=2, xpd=TRUE)
    text(Rax[4,1],Rax[4,2], labels="y", pos=4, xpd=TRUE)
    text(Rax[6,1],Rax[6,2], labels="z", pos=3, xpd=TRUE)
    

    segments(Rbox[1,1],Rbox[1,2], Rbox[2,1],Rbox[2,2], col=boxcol, xpd=TRUE)
    segments(Rbox[2,1],Rbox[2,2], Rbox[3,1],Rbox[3,2], col=boxcol, xpd=TRUE)
    segments(Rbox[3,1],Rbox[3,2], Rbox[4,1],Rbox[4,2], col=boxcol, xpd=TRUE)
    segments(Rbox[4,1],Rbox[4,2], Rbox[1,1],Rbox[1,2], col=boxcol, xpd=TRUE)

    segments(Rbox[5,1],Rbox[5,2], Rbox[6,1],Rbox[6,2], col=boxcol, xpd=TRUE)
    segments(Rbox[6,1],Rbox[6,2], Rbox[7,1],Rbox[7,2], col=boxcol, xpd=TRUE)
    segments(Rbox[7,1],Rbox[7,2], Rbox[8,1],Rbox[8,2], col=boxcol, xpd=TRUE)
    segments(Rbox[8,1],Rbox[8,2], Rbox[5,1],Rbox[5,2], col=boxcol, xpd=TRUE)

    segments(Rbox[2,1],Rbox[2,2], Rbox[6,1],Rbox[6,2], col=boxcol, xpd=TRUE)
    segments(Rbox[3,1],Rbox[3,2], Rbox[7,1],Rbox[7,2], col=boxcol, xpd=TRUE)
    segments(Rbox[4,1],Rbox[4,2], Rbox[8,1],Rbox[8,2], col=boxcol, xpd=TRUE)
    segments(Rbox[1,1],Rbox[1,2], Rbox[5,1],Rbox[5,2], col=boxcol, xpd=TRUE)
  }

