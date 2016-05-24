getseis24<-function(DB, iyear=2009, iday=1, usta="", acomp="", kind = 1,  Iendian=1, BIGLONG=FALSE )
  {
    if (missing(kind)) {
      kind = 1
    }
      if(missing(Iendian)) { Iendian=1 }
     if(missing(BIGLONG)) { BIGLONG=FALSE}

    origyr = min(DB$yr, na.rm = TRUE)
    attr(DB, "origyr")<- origyr
    eday = EPOCHday(iyear, jd = iday, origyr = origyr)

    M1 = eday$jday

    at1 = M1 +0/24
    at2 = at1 + 24/24
    


    w1 = which(at1 >= DB$t1 & at1 < DB$t2)
    w2 = which(at2 >= DB$t1 & at2 < DB$t2)
    w3 = which(DB$t1 >= at1 & DB$t1 < at2)
    w4 = which(DB$t2 >= at1 & DB$t2 < at2)

    if (length(c(w1, w2, w3, w4)) < 1) {
      print("getseis24: No Match in DataBase")
      return()
    }
    wi = unique(c(w1, w2, w3, w4))
    
##     wi = unique(c(w1, w2, w3, w4))
    if (length(wi) < 1) {
      print("getseis24: No times match this call")
      return()
    }
    fnloc = DB$fn[wi]

    fn1 = paste(sep=".", DB$fn[wi], DB$sta[wi], DB$comp[wi])


    
    
    nsta = length(usta)
    if (is.null(acomp))
      acomp = "*"
    ncomp = length(acomp)
    gi = vector()
    for (i in 1:nsta) {
      for (j in 1:ncomp) {
        gi = c(gi, grep(paste(sep = ".", usta[i], acomp[j]),
          fn1))
      }
    }
    if (length(gi) < 1) {

      print("getseis24: length(gi) < 1")
      return()
    }

    
    fn2 = fnloc[gi]

    ###   read in the data from the SEGY data base and store
    GG = JGET.seis(fn2, kind = kind, PLOT = -1, Iendian=Iendian, BIGLONG=BIGLONG )
###  this set of data should have only one station and one component

    ###  set up some information from each of the traces in the data set
    N = length(GG)
    stas = rep(NA, N)
    comp = rep(NA, N)
    units = rep(NA, N)
    KNOTE = rep(NA, N)
    gdt = rep(NA, N)
    gt1 = rep(NA, N)
    gt2 = rep(NA, N)
    gdur = rep(NA, N)
    LENS = rep(NA, N)
    YRS = rep(NA, N)
    for (i in 1:N) {
      stas[i] = GG[[i]]$sta
      comp[i] = GG[[i]]$comp
      units[i] = GG[[i]]$units
      KNOTE[i] = paste(sep = ".", GG[[i]]$sta, GG[[i]]$comp)
      gdt[i] = GG[[i]]$dt
      YRS[i] = GG[[i]]$DATTIM$yr
      eday1 = EPOCHday(GG[[i]]$DATTIM$yr, jd = GG[[i]]$DATTIM$jd, origyr = origyr)
      gt1[i] = eday1$jday + GG[[i]]$DATTIM$hr/24 + GG[[i]]$DATTIM$mi/(24 *
          60) + GG[[i]]$DATTIM$sec/(24 * 3600) + GG[[i]]$DATTIM$msec/(1000 *
                                                                      24 * 3600)
      gdur[i] = (GG[[i]]$N * GG[[i]]$dt)/(24 * 3600)
      LENS[i] = GG[[i]]$N
    }

    #################

    gt2 = gt1 + gdur
    
    adt = gdt[1]
    FIX = 24
    h = FIX


   ##    set up structure for output
    GY = list(
      yr=rep(iyear, length=24),
      jd = rep(iday, length=24),
      t1=M1+ (0:23)/24,
      t2=M1+ (1:24)/24,
      ed=rep(M1, length=24),
      hr=   0:23,
      mi=  rep(0, length=24), 
      sec=rep(0, length=24),
      gamp=rep(0, length=24),
      gdt=rep(adt, length=24),
      gnam=rep(KNOTE[1], length=24),
      gfile=rep(iyear, length=24) ,
      sigs=as.list(1:24),
      zna = as.list(1:24)  )
    

  
    
    xa = seq(from=0, length=3600/adt, by=adt)
    
    for(i in 1:h )
      {
        N2 = (i%%2)+1
###  here make a minor adjustment to avoid sampling problems
        a1 = M1 + (i-1)/24 - 0*adt
        a2 = M1 + (i)/24 +  0*adt
        
        ##print(paste(a1, a2))
        w1 = which(gt1>=a1&gt1<=a2)
        w2 = which(gt2>=a1&gt2<=a2)
        w = sort(unique(c(w1,w2)))
        
        if(length(w)>1)
          {
            zed = rep(NA, length=3600/adt)
            ex =  a1+xa/(24*3600)
            
            for(j in 1:length(w))
              {
                k = w[j]
                t1 = gt1[k]
                t2 = gt2[k]
                f1 = ex>=t1&ex<=t2
###    print(paste(sep=' ', j, length(f1) )) 
                U = any(f1)
                if(U)
                  {
                    s = GG[[k]]$amp
                    if(gdt[k]==adt)
                      {
                        
                        lex = gt1[k]+seq(from=0, by=gdt[k], length=(length(s)))/(24*3600)
                        tem = lex>=a1&lex<=a2
                        ## print(paste(sep=' ', j, length(zed[f1]) , length(s[tem])))
                        s2 = s[tem]
                        zed[f1] = s2[1:length(zed[f1])]
                      }
                    else
                      {  ###########  need to  resample, I guess
                        
                        
                        lex = gt1[k]+seq(from=0, by=gdt[k], length=(length(s)))/(24*3600)
                        
                        
                        slex = spline(lex, s, n = length(xa) , method = "fmm",
                          xmin = min(lex), xmax = max(lex), ties = mean)
                        
                        
                        tem = slex$x>=a1&slex$x<=a2
                        ## print(paste(sep=' ', j, length(zed[f1]) , length(s[tem])))
                        s2 = slex$y[tem]
                        zed[f1] = s2[1:length(zed[f1])]
                        

                      }
                    ## print(U)
                  }
              }




            
            y1 = -a1
            amean = mean(zed, na.rm=TRUE )
            
            ##   zed=zed-amean
            
            
            zna = is.na(zed)
            
            if(any(zna))
              {
                zed[is.na(zed)] = amean
              }
            
            
            GY$sigs[[i]] = zed
            GY$zna[[i]]  = zna   
            GY$yr[i] = iyear
            GY$t1[i]=a1
            GY$t2[i]=a2
            GY$ed[i]=M1
            GY$hr[i]=i-1
            GY$sec[i]=0
            GY$dt[i]=adt
            
            
            
            
            
            
            
          }


        
      }

    
    invisible(GY)
     
        
        

  }
