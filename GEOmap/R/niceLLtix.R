niceLLtix<-function(rcoords)
  {
    
#####  examples
   
###   rcoords = c(91.5, 93.8)

    K = 10

    scoords = sort(range(rcoords) )


    dms1 = dms(scoords)
    
    if(diff(dms1$d)<=1)
      {
        ######  same degree, check minutes
        dm = diff(dms1$m)
        
        if(abs(dm)<=1)
          {


            dsec = diff(dms1$m*60+dms1$s)
             K = goodticdivs(dsec)
            secstart = K*trunc(dms1$s[1]/K)
            secend = secstart+60*(dm+1)
            print(paste(sep=" ","sec", secstart, secend, K))
                  
            isec = seq(from=secstart, to=secend, by=K)
            
            FDEGcoord = rep(dms1$d[1], times=length(isec))
            FMINcoord = rep(dms1$m[1], times=length(isec))
            FSECcoord = isec


              
          }
        else
          {

            dmin = diff(dms1$d*60+dms1$m)

              K = goodticdivs(dmin)
          
             
            
            minstart = K*trunc(dms1$m[1]/K)
            minend = minstart+(dmin+1)

            
             print(paste(sep=" ", "min", minstart, minend, K))
            imin = seq(from=minstart, to=minend, by=K)
            
            FDEGcoord = rep(dms1$d[1], times=length(imin))
            FMINcoord = imin
            FSECcoord =  rep(0, times=length(imin))

            }

          }
    else
      {
        ######  degrees vary considerable

        
        ddeg = diff(dms1$d)

        K = 15
        
        K = goodticdivs(ddeg)

            
            degstart = K*trunc(dms1$d[1]/K)
            degend =degstart+(ddeg+1)

              print(paste(sep=" ", "deg", degstart, degend, K))
            ideg = seq(from=degstart, to=degend, by=K)
            
            FDEGcoord =  ideg
            FMINcoord =  rep(0, times=length(ideg))
            FSECcoord =  rep(0, times=length(ideg))

        
      }
       
     
    DDcoord = FDEGcoord+FMINcoord/60+FSECcoord/3600
    w = which(DDcoord>=scoords[1] & DDcoord<=scoords[2])
    ## FDEGcoord =  FDEGcoord[w]
   ##  FMINcoord =  FMINcoord[w] 
    ## FSECcoord =  FSECcoord[w]
    DDcoord =DDcoord[w]


    si = sign(DDcoord)

    final = dms(abs(DDcoord))
    
     ##   return(list(DD=DDcoord, deg=FDEGcoord, min=FMINcoord, sec=FSECcoord))
  

    return(list(DD=DDcoord, deg=final$d, min=final$m, sec=final$s, si=si))
  
  }

