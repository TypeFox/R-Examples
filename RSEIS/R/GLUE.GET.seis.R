`GLUE.GET.seis` <-
function(GG)
      {
          ####  names(GG[[1]])
     ####   print("working in GLUE.GET.seis")
        
	N = length(GG)
        stas  = rep(NA, N)
        comp =  rep(NA, N)
        units =  rep(NA, N)
        
        KNOTE = rep(NA, N)
        dt =  rep(NA, N)
        t1 =  rep(NA, N)
        t2 =  rep(NA, N)
        dur =  rep(NA, N)
        
        LENS =  rep(NA, N)
        YRS =  rep(NA, N)

        
        

        ##############  find the minimum start time for all traces
        #######   if traces cross a year boundary need to make special
        #######   case

        ## gather information here

        ## t1 and t2 are in julian days - but if the year boundary is crossed,
        ##  this will be wrong

        for(i in 1:N)
          {
            stas[i] = GG[[i]]$sta
            comp[i] = GG[[i]]$comp
            units[i] = GG[[i]]$units
            
            KNOTE[i] = paste(sep=".",GG[[i]]$sta, GG[[i]]$comp)
            dt[i] = GG[[i]]$dt
            YRS[i] =  GG[[i]]$DATTIM$yr


            eday = EPOCHday(GG[[i]]$DATTIM$yr, jd=GG[[i]]$DATTIM$jd)
            
            t1[i] = eday$jday+GG[[i]]$DATTIM$hr/24+GG[[i]]$DATTIM$mi/(24*60)+
              GG[[i]]$DATTIM$sec/(24*3600)+GG[[i]]$DATTIM$msec/(1000*24*3600)

             LENS[i] = length(GG[[i]]$amp)
            dur[i] =  (LENS[i]*GG[[i]]$dt)/(24*3600)

            
            
           
          }



 
        t2 = t1 + dur
        
        ############  combine traces that have the same station name and component

        UK = unique(KNOTE)

        RR = vector(mode="list")

        for(i in 1:length(UK))
          {
            j = which(match( KNOTE, UK[i])==1)

            w1 = which.min(t1[j])
            
            w2 = which.max(t2[j])

            mint = t1[j[w1]]

            ###  this combines the length of all matches

          ### print(c(i, j[w1], t1[j[w1]], t2[j[w1]], j[w2],  t1[j[w2]], t2[j[w2]] ))

                  
            temy = rep(NA, length(seq(from=0, to=(t2[j[w2]]-mint)*24*3600  , by=dt[j[w1]])))
            
            for( k in j)
              {
                samp1 = floor( ((t1[k]-mint)*24*3600)/dt[k])+1
                
                samp2 =  samp1+LENS[k]-1
                
                temy[samp1:samp2] = GG[[k]]$amp

              }

            nsamp=length(temy)
            tem2 = nsamp*dt[j[w1]]

            if(is.null(GG[[j[w1]]]$HEAD$scale_fac))
              { scalefac = 1 } else {scalefac =GG[[j[w1]]]$HEAD$scale_fac} 
            
            if(is.null(GG[[j[w1]]]$HEAD$gainConst))
              { gain = 1 } else {gain = GG[[j[w1]]]$HEAD$gainConst} 

       

            
            RR[[i]] = list(fn=GG[[j[w1]]]$fn, sta=GG[[j[w1]]]$sta,
                units=GG[[j[w1]]]$units , comp=GG[[j[w1]]]$comp, dt=GG[[j[w1]]]$dt,
                DATTIM=list(yr=GG[[j[w1]]]$DATTIM$yr, jd=GG[[j[w1]]]$DATTIM$jd,
                  mo=GG[[j[w1]]]$DATTIM$mo,  dom=GG[[j[w1]]]$DATTIM$dom,
                  hr=GG[[j[w1]]]$DATTIM$hr ,   mi=GG[[j[w1]]]$DATTIM$mi ,
                  sec=GG[[j[w1]]]$DATTIM$sec ,  msec=GG[[j[w1]]]$DATTIM$msec ,
                  dt=GG[[j[w1]]]$DATTIM$dt ,
                  t1=GG[[j[w1]]]$DATTIM$t1 ,  t2=tem2 ,   off=GG[[j[w1]]]$DATTIM$off),
                N=nsamp ,
                gain=gain, scalefac=scalefac, 
                amp=temy)
            
          }


        invisible(RR)
      }

