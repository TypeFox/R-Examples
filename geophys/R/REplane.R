REplane <-
function(m1, Lp, PPs, Rbox, Rview, xscale)
  {

  
    if(m1==1)
      {
        points(Lp)
        basepoint = 3
##########   point 1: check for line closest to either 7-3, 4-3, 2-3

        legpoints = c(7,4,2)

     ###   annotatebox()
   ###      Lp = locator(1)
        ###  m1 =1

       
        

        VL =   cbind( rep(Rbox[basepoint,1] , length(legpoints)), rep(Rbox[basepoint,2] , length(legpoints)),
          Rbox[legpoints,1], Rbox[legpoints,2])

        G = points2line(Lp, VL )

        w1 = which.min(G$srat)

        rat = G$rat[w1]
        
       ##
        if(w1==1)
          {
            PPs[m1, 1]  =  xscale*rat
            PPs[m1, 2]  =  xscale*1
            PPs[m1, 3]  =  xscale*1           
          }

        if(w1==2)
          {
            PPs[m1, 1]  =  0
            PPs[m1, 2]  =  xscale*(1-rat)
            PPs[m1, 3]  =  xscale*1          
          }

        if(w1==3)
          {
            PPs[m1, 1]  =  0
            PPs[m1, 2]  =  xscale*1
            PPs[m1, 3]  =  xscale*(1-rat)         
          }




        
        
      }
    if(m1==2)
      { 
        
        points(Lp)
        basepoint = 8
##########   point 1: check for line closest to either 7-3, 4-3, 2-3

        legpoints = c(7, 4, 5)

  
        VL =   cbind( rep(Rbox[basepoint,1] , length(legpoints)), rep(Rbox[basepoint,2] , length(legpoints)),
          Rbox[legpoints,1], Rbox[legpoints,2])

        G = points2line(Lp, VL )

        w1 = which.min(G$srat)
        rat = G$rat[w1]
        
       ##
        if(w1==1)
          {
            PPs[m1, 1]  =  xscale*1 
            PPs[m1, 2]  =  xscale*rat
            PPs[m1, 3]  =  xscale*1           
          }

        if(w1==2)
          {
            PPs[m1, 1]  =  xscale*(1-rat)
            PPs[m1, 2]  =  0
            PPs[m1, 3]  =  xscale*1          
          }

        if(w1==3)
          {
            PPs[m1, 1]  =  xscale*1
            PPs[m1, 2]  =  0
            PPs[m1, 3]  =  xscale*(1-rat)         
          }



      }
    if(m1==3)
      {
        points(Lp)
        basepoint = 6
##########   point 1: check for line closest to either 7-3, 4-3, 2-3

        legpoints = c(7, 5, 2)

  
        VL =   cbind( rep(Rbox[basepoint,1] , length(legpoints)), rep(Rbox[basepoint,2] , length(legpoints)),
          Rbox[legpoints,1], Rbox[legpoints,2])

        G = points2line(Lp, VL )

        w1 = which.min(G$srat)
        rat = G$rat[w1]
        
       ##
        if(w1==1)
          {
            PPs[m1, 1]  =  xscale*1 
            PPs[m1, 2]  =  xscale*1
            PPs[m1, 3]  =  xscale*rat        
          }

        if(w1==2)
          {
            PPs[m1, 1]  =  xscale*1
            PPs[m1, 2]  =   xscale*(1-rat)
            PPs[m1, 3]  =  0         
          }

        if(w1==3)
          {
            PPs[m1, 1]  =  xscale*(1-rat)
            PPs[m1, 2]  =  xscale*1
            PPs[m1, 3]  =  0         
          }


      }
    return(PPs)
  }

