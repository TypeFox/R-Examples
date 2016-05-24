readCMT <-
function(filename, PLOT=TRUE)
  {

############  read and parse a file from the Harvard CMT database
#########   in  CMTSOLUTION format

###########  data seems to separated by a blank line between records
    if(missing(PLOT)) { PLOT=TRUE }

    
    read1CMTrecord<-function(con1)
      {
        OneCMT = vector( mode='character', length = 0 )
        i = 0
        r1  = 0
        while( !identical(r1, "") )
          {
            i = i+ 1
            r1 = readLines(con1, n=1)
            if(is.null(r1)) { break }
            if(length(r1)<1) { break }
            
            OneCMT[i] = r1

          }
        
        return(OneCMT)
        
      }

    ParseCMTrecord<-function(H1)
      {
        Hlist = list()
        
        L1 = unlist( strsplit(H1[1], split=" "))
        L1 = L1[L1!=""]

        Hlist$FIRST = H1[1]


        ####### the files are screwed up - the year is sometimes attached to PDEW
        #######   needs care in handling

        if( as.numeric(L1[2])> 12)
          {
            Hlist$yr = as.numeric(L1[2])
            Hlist$mo = as.numeric(L1[3])
            Hlist$dom = as.numeric(L1[4])
            Hlist$hr  = as.numeric(L1[5])
            Hlist$mi  = as.numeric(L1[6])
            Hlist$sec  = as.numeric(L1[7])
            

          }
        else
          {
            tem = L1[1]
             Hlist$yr = as.numeric(substr(tem, 5, 8))
            Hlist$mo = as.numeric(L1[2])
            Hlist$dom = as.numeric(L1[3])
           
            Hlist$hr  = as.numeric(L1[4])
            Hlist$mi  = as.numeric(L1[5])
            Hlist$sec  = as.numeric(L1[6])

          }

        i1 = grep("event name:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$name = L2[3]
        i1 = grep("time shift:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$tshift = as.numeric(L2[3])
        i1 = grep("half duration:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$half = L2[3]
        i1 = grep("latitude:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$lat = as.numeric(L2[2])
        i1 = grep("longitude:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$lon = as.numeric(L2[2])
        i1 = grep("depth:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$z = as.numeric(L2[2])
        i1 = grep("Mrr:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$Mrr = as.numeric(L2[2])
        i1 = grep("Mtt:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$Mtt = as.numeric(L2[2])
        i1 = grep("Mpp:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$Mpp = as.numeric(L2[2])
        i1 = grep("Mrt:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$Mrt = as.numeric(L2[2])
        i1 = grep("Mrp:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$Mrp = as.numeric(L2[2])
        i1 = grep("Mtp:", H1)
        L2   =  unlist( strsplit(H1[i1], split=" "))
        L2 = L2[L2!=""]
        Hlist$Mtp = as.numeric(L2[2])
        
        return(Hlist)
      }

    if(missing(filename)) filename = "CMT_FULL_FORMAT.txt"

    
    con1 = file(filename, open = "r")

    k = 0
    H1 = c(1,2, 3)
    P1 = NULL

    GetCMT = list()
    
    while(length(H1)>2)
      {
        H1 = read1CMTrecord(con1)
        if(length(H1)>6) P1 = ParseCMTrecord(H1)
        
###### print(length(H1))
######  print(H1)
       
        if(length(P1$Mrr)>0)
          {
            k = k + 1
            GetCMT[[k]] = P1

            if(PLOT)
              {
                moments = cbind(k, P1$Mtt, P1$Mpp, P1$Mrr, -P1$Mrp, P1$Mrt ,-P1$Mtp)

                jo = doNonDouble(moments, sel=1)
                title(sub=P1$name)
                
                locator()
              }
          }

      }


    
    
    close(con1)
    
    
    invisible(GetCMT)

  }
