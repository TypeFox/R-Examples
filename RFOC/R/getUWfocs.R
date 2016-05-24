`getUWfocs` <-
function(amfile)
  {
##############   FS = getAMcards(amfile)
   ####  require(RSEIS)

    getMCARD<-function(MCARD)
      {
           F = list(az=0, dip=0)
          F$az = as.numeric(substr(MCARD, 5, 7))
          F$dip = as.numeric(substr(MCARD,8, 10))
          
          G = list(az=0, dip=0)
          G$az = as.numeric(substr(MCARD, 14, 16))
          G$dip = as.numeric(substr(MCARD,17, 19))
          
          U = list(az=0, dip=0)
          U$az = as.numeric(substr(MCARD, 22, 25))
          U$dip = as.numeric(substr(MCARD, 26, 28))
          
          V = list(az=0, dip=0)
          V$az = as.numeric(substr(MCARD, 32, 34))
          V$dip = as.numeric(substr(MCARD, 35, 37))
          
          P = list(az=0, dip=0)
          P$az = as.numeric(substr(MCARD, 41, 43))
          P$dip = as.numeric(substr(MCARD, 44, 46))
          
          T  = list(az=0, dip=0)
          T$az = as.numeric(substr(MCARD, 50, 52))
          T$dip = as.numeric(substr(MCARD, 53, 55))
          
          
          MC = list(F=F, G=G, U=U, V=V, P=P, T=T, CNVRG=NA)

           return(MC)
        
        
      }
    
    


    CMT = list(lon=0, lat=0, str1=0, dip1=0, rake1=0, str2=0, dip2=0, rake2=0, sc=0, iexp=0, name="", yr=0, mo=0, dom=0, jd=0, hr=0, mi=0, se=0)

 
    
    zz <- file(amfile, "r")
    ##  open(zz)
    APF = readLines(con = zz, ok = TRUE, warn = TRUE)
    close(zz)
    N = length(APF)

    K = list()
    j = 0
    for(i in seq(from=1, to=N, by=2))
      {
       
        
        GA =  RSEIS::unpackAcard(APF[i])
        GF =  getMCARD(APF[i+1])
        
        K[[i]] = list(A=GA, M=GF)
        j = j+1
        CMT$yr[j] = GA$yr
        CMT$mo[j] = GA$mo
        CMT$dom[j] = GA$dom
        CMT$jd[j] =  GA$jd
        CMT$hr[j] = GA$hr
        CMT$mi[j] = GA$mi
        CMT$se[j] = GA$se
        
        CMT$lat[j] = GA$lat
        CMT$lon[j] = GA$lon
        CMT$z[j] = GA$z
        
        CMT$mag[j] = GA$mag
        
        
        
        ang2 = GetRakeSense(GF$U$az, GF$U$dip, GF$V$az, GF$V$dip, GF$P$az, GF$P$dip,  GF$T$az, GF$T$dip)

        MEC = GetRake(GF$F$az-90, GF$F$dip,  GF$G$az -90,  GF$G$dip, ang2)


        CMT$str1[j] = GF$F$az-90
        CMT$dip1[j] = GF$F$dip
        CMT$rake1[j] = MEC$rake1

        CMT$str2[j] = GF$G$az-90
        CMT$dip2[j] = GF$G$dip
        CMT$rake2[j] = MEC$rake2
        
        
      }
   
 invisible(CMT)

   
  }

