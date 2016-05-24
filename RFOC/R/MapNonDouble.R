MapNonDouble <-
  function(Locs, moments, sel=1, siz = 0.2, col=rgb(1, .75, .75) , PLANES = TRUE, add=FALSE, LEG=FALSE )
  {
###   originally developed for matlab by Vaclav Vavrycuk
### see: http://www.ig.cas.cz/en/research-&-teaching/software-download/
### translated to R by Jonathan M. Lees
    if(missing(col)) { col =rgb(1, .75, .75) }

    if(missing(PLANES)) {  PLANES = TRUE }
    if(missing(LEG)) { LEG=FALSE  }
  
    
    di = dim(moments)
    number.of.events = di[1]
    moment_11 = moments[,2]
    moment_22 = moments[,3]
    moment_33 = moments[,4]
    moment_23 = moments[,5]
    moment_13 = moments[,6]
    moment_12 = moments[,7]

    if(missing(sel)) { sel = 1:length(Locs$x) }
    if(missing(add)) {add=FALSE }
    if(missing(siz)) {    siz = 0.2 }  ###  size of glyph in inches

    

    if(!add) plot(range(Locs$x), range(Locs$y), asp=1, xlab="EW", ylab="NS" )

                                        # points(Locs$x,Locs$y)

    up = par("usr")
    ui = par("pin")

    ratx = (up[2]-up[1])/ui[1]
    raty=  (up[4]-up[3])/ui[2]

    
    
###graphics.off()
###dev.new()
    

    angles.all = matrix(0, ncol=6, nrow=length(sel))
    LAB.all = matrix(0, ncol=2, nrow=length(sel))

    Fi=seq(from=0, by=0.1, to=361)

    CX = cos(Fi*pi/180.0)
    CY = sin(Fi*pi/180.0)
    
### ------------------------------------------------
###  the main loop over events
    for (j  in 1:length(sel)  )
      {
        i = sel[j]
        m=matrix( c(moment_11[i],moment_12[i],moment_13[i],
          moment_12[i],moment_22[i],moment_23[i],
          moment_13[i],moment_23[i],moment_33[i]), ncol=3, byrow=TRUE)
###    print(m)
### -----------------------------------------
###  calculation of a fault normal and slip from the moment tensor
        angles.all[j,] = FOCangles(m)
        strike.1 = angles.all[j,1]
        dip.1    = angles.all[j,2]
        rake.1   = angles.all[j,3]
        
### -------------------------------------
###  plotting the solution
        
        usizx = siz*ratx
        usizy = siz*raty


        polygon(Locs$x[i]+usizx*CX,Locs$y[i]+usizy*CY, col='white')

        if(is.null(col))
          { jcol = foc.color(foc.icolor(rake.1), pal=1) }
        else
          {
            jcol = col

          }

        
        shad =  ShadowCLVD(m, PLOT=FALSE)

        if(!all(is.na(shad$z)))
          {
        image(Locs$x[i]+usizx*shad$x,   Locs$y[i]+usizy*shad$y, z=shad$z, col=jcol , add=TRUE )
      }
        
        lines(Locs$x[i]+usizx*CX,Locs$y[i]+usizy*CY)
        LAB.all[j,] = c(Locs$x[i]+usizx*CX[1], Locs$y[i]+usizy*CY[1])
        
        if(PLANES)
          {
            NODS = nodalLines(strike.1, dip.1, rake.1, PLOT=FALSE)

            lines(Locs$x[i]+usizx*NODS$PLANE1$x,   Locs$y[i]+usizy*NODS$PLANE1$y)
            lines(Locs$x[i]+usizx*NODS$PLANE2$x,   Locs$y[i]+usizy*NODS$PLANE2$y)
          }

    
      }


    if(LEG)
      {
        
        fleg = 1:7
        flegc = foc.color(fleg, pal=1)
        flab = focleg(fleg)

        legend("topleft", legend=flab, col=flegc, pch=19, bg="white" )
      }



invisible(list(FOC=angles.all, LAB=LAB.all) ) 
}
