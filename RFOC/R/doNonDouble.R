doNonDouble <-
function(moments, sel=1, col=rgb(1, .75, .75))
  {
###   originally developed for matlab by Vaclav Vavrycuk
### see: http://www.ig.cas.cz/en/research-&-teaching/software-download/
### translated to R by Jonathan M. Lees
    if(missing(col)) { col =rgb(1, .75, .75) }
    
    di = dim(moments)
    number.of.events = di[1]
    moment_11 = moments[,2]
    moment_22 = moments[,3]
    moment_33 = moments[,4]
    moment_23 = moments[,5]
    moment_13 = moments[,6]
    moment_12 = moments[,7]


    if(missing(sel)) { sel = 1:number.of.events }


    angles.all = matrix(0, ncol=6, nrow=length(sel))

    
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
### ----------------------------------------
###  boundary RPMG::circle
        Fi=seq(from=0, by=0.1, to=361)
###  dev.new()
        plot(cos(Fi*pi/180.0),sin(Fi*pi/180.0),type='l', asp=1 , ann=FALSE, axes=FALSE)
        
        ShadowCLVD(m, col=col)
###      shadowing(m)


###  denoting the North direction
        segments( 0, 1,  0, 1.07 )
        text(0.0, 1.05,'N', pos=3, xpd=TRUE)
        
        Nod = nodalLines(strike.1, dip.1, rake.1, PLOT=TRUE)

        
###### strike= strike.1; dip=dip.1; rake=rake.1
        PTaxes(strike.1, dip.1, rake.1)
        title(main=paste("event=", i))
        if(j<length(sel)   ) locator(1)
      }

invisible( angles.all) 
}
