stressSETUP <-
function(P1=c(.2, 1,1,0 ), P2=c(1, .1,1,0), P3=c(1, 1,.4,0), xscale=30)
  {


    if(missing(xscale)) { xscale=30 }
    if(missing(P1)) P1=c(.2, 1,1,0)
    if(missing(P2)) P2=c(1, .1,1,0)
    if(missing(P3)) P3=c(1, 1,.4,0)

    if(length(P1)<4) P1 = c(P1, 0)
    if(length(P2)<4) P2 = c(P2, 0)
    if(length(P3)<4) P3 = c(P3, 0)
    

    #### this is basic setup for the box rotation routines
    
   ####  require(RFOC)
    
    Rview  =    RFOC::ROTZ(-125) %*% RFOC::ROTX(-55) 

    BOX <-matrix(c(0,0,0,0,
                   0, 1, 0,0,
                   0, 1, 1,0,
                   0, 0, 1,0,
                   1,0,0,0,
                   1, 1, 0,0,
                   1, 1, 1,0,
                   1, 0, 1,0), ncol=4, byrow=TRUE)



    BOX = xscale*BOX


    AX = matrix(c(0,0,0,0,
      1, 0, 0,0,
      0, 0, 0,0,
      0, 1, 0,0,
      0,0,0,0,
      0, 0, 1,0), ncol=4, byrow=TRUE)

    AX = 1.5*xscale*AX


    Rax =  AX %*% Rview


    Rbox =   BOX %*% Rview

    headlen =xscale* .3/6
    len =xscale* .7/6
    basethick =xscale* 0.05/2
    headlip =xscale* .02/2
    aglyph = RFOC::Z3Darrow(len = len , basethick =basethick , headlen =headlen , headlip=headlip )

    
    P1 = xscale*P1
    P2 = xscale*P2
    P3 = xscale*P3

    PPs = rbind(P1, P2, P3)
    
    Rp = PPs  %*% Rview


    invisible(list(xscale=xscale, Rview=Rview, BOX=BOX, AX=AX, Rbox=Rbox, Rax=Rax, PPs=PPs, Rp=Rp, aglyph=aglyph))
}

