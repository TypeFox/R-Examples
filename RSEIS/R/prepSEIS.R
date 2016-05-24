`prepSEIS` <-
function(GG)
{
  ###  prepare a list of seismic information after reading in segy/sac/ah files
  ###    using SEE.ahseis(fnames, kind=1, PLOT=TRUE)
  gstas = rep(NA, length(GG))
  gcomps = rep(NA, length(GG))
  
  gtim1 = rep(NA, length(GG))
  gtim2 = rep(NA, length(GG))
  gn = rep(NA, length(GG))
  gdt  = rep(NA, length(GG))
    ## gfn  = rep(NA, length(GG))


  ####  prepare some of the stats on the times of the waveforms
  oyears = rep(1972, length(GG))
  
  for(i in 1:length(GG)) {  oyears[i] = GG[[i]]$DATTIM$yr  }
  eyear  =min(oyears, na.rm=TRUE)
  
  for(i in 1:length(GG))
    {
      ### gstas[i] = GG[[i]]$sta
      dt = round(100000*GG[[i]]$DATTIM$dt)/100000
      n = length(GG[[i]]$amp)

      if(is.null(GG[[i]]$DATTIM$jd))
         {
           jd =getjul(GG[[i]]$DATTIM$yr, GG[[i]]$DATTIM$mo,GG[[i]]$DATTIM$dom )
           GG[[i]]$DATTIM$jd = jd

         }
      if(is.null(GG[[i]]$DATTIM$mo))
         {
           md = getmoday(GG[[i]]$DATTIM$jd, GG[[i]]$DATTIM$yr)

           GG[[i]]$DATTIM$mo = md[1]
           GG[[i]]$DATTIM$dom = md[2]
                    
         }

      if(is.null(GG[[i]]$gain))
         {
           GG[[i]]$gain = 1
         }
     if(is.null(GG[[i]]$scalefac))
         {
           GG[[i]]$scalefac = 1
         }

      GG[[i]]$DATTIM$dt = dt

      eday = EPOCHday(GG[[i]]$DATTIM$yr, jd=GG[[i]]$DATTIM$jd, origyr = eyear )
      
      gtim1[i]  = eday$jday+GG[[i]]$DATTIM$hr/24+GG[[i]]$DATTIM$mi/(24*60)+GG[[i]]$DATTIM$sec/(24*3600)
      gtim2[i]  = gtim1[i]+GG[[i]]$N*GG[[i]]$DATTIM$dt/(24*3600)
      gn[i] = n
      gdt[i] = dt
      ## gfn[i] = GG[[i]]$fn
   ##  print(paste(sep=' ', i, GG[[i]]$sta, GG[[i]]$comp, GG[[i]]$N, n, GG[[i]]$DATTIM$dt, GG[[i]]$DATTIM$jd, GG[[i]]$DATTIM$hr, GG[[i]]$DATTIM$mi, GG[[i]]$DATTIM$sec))
     ### print(paste(sep=' ', i, dt, GG[[i]]$DATTIM$dt))
    }

  if(any(is.infinite(gdt))) return(NULL)
  
    wmin = which.min(gtim1)
    wmax = which.max(gtim2)

####  set up the padding
  r1 =  round((gtim1-gtim1[wmin])*24*3600/gdt)
  r2 =  round((gtim2[wmax]-gtim2)*24*3600/gdt)

  r2[r2<1] = 0
  r1[r1<1] = 0
  r2[is.nan(r2)] = 0
  r2[is.infinite(r2)] = 0
  
  BigR = r1+gn+r2
  ma = 1:length(gdt)
   
  K = BigR[ma[1]]
  ascd = as.list(1:length(GG))
   
  dt = gdt
  notes = rep(NA, length(ma))
  stns = rep(NA, length(ma))
  comps = rep(NA, length(ma))
  units = rep(NA, length(ma))
  
  info = list(fn=rep(NA, length(ma)), name=rep(NA, length(ma)), yr=rep(NA, length(ma)), jd=rep(NA, length(ma)),
    mo=rep(NA, length(ma)), dom=rep(NA, length(ma)), hr=rep(NA, length(ma)), mi=rep(NA, length(ma)),
    sec=rep(NA, length(ma)), msec=rep(0, length(ma)),dt=rep(0, length(ma)),t1=rep(0, length(ma)) ,t2=rep(0, length(ma)),
    off=rep(0, length(ma)), n1=rep(0, length(ma)), n2=rep(0, length(ma)), n3=rep(0, length(ma)), n=rep(0, length(ma)),
    gain=rep(1, length(ma))   , scalefac=rep(1, length(ma))
 
    )
  ###  fill up data structure and information
  ############   the NA padding of traces occurs here
  ###   source("/home/lees/Progs/R_stuff/seis.R"); save.image()

  for(j in 1:length(ma))
	{
	  ima = ma[j]
     ###     print(paste(sep=" ", j, ima, r1[ ima ], r2[ ima ]))
          ######   the NA padding of traces occurs here
          
	  ascd[[j]] = c(rep(NA, r1[ ima ]) , GG[[ima]]$amp,    rep(NA,r2[ ima ]))

	  notes[j] = paste(sep=' ', GG[[ima]]$sta, GG[[ima]]$comp)
	  stns[j] = GG[[ima]]$sta
	  comps[j] = GG[[ima]]$comp
          units[j] = "volts"
          
        ##   print(paste(sep=' ', GG[[ima]]$sta, GG[[ima]]$comp, GG[[ima]]$units))
	 
	  info$fn[j] = GG[[ima]]$fn
	  info$name[j] = GG[[ima]]$fn
	  info$yr[j] = GG[[ima]]$DATTIM$yr
	  info$jd[j] = GG[[ima]]$DATTIM$jd
	  info$mo[j] = GG[[ima]]$DATTIM$mo
	  info$dom[j] = GG[[ima]]$DATTIM$dom
	  info$hr[j] = GG[[ima]]$DATTIM$hr
	  info$mi[j] = GG[[ima]]$DATTIM$mi
	  info$sec[j] = GG[[ima]]$DATTIM$sec
	  info$msec[j] = 0
	  info$dt[j] = gdt[ima]

 	  info$t1[j] = 0
 	  info$t2[j] = gdt[ima]*(length(ascd[[j]])-1)
           
	  info$off[j] =  r1[ ima ]*gdt[ima]
	  info$n1[j] =  length(ascd[[j]])
	  info$n2[j] =  info$n1[j]
 	  info$n3[j] =  info$n1[j]
 	  info$n[j] =  info$n1[j]

 	  info$gain[j] =  GG[[ima]]$gain

 	  info$scalefac[j] = GG[[ima]]$scalefac
          
	}
 
    f1 = unlist(strsplit(info$fn[1], "/"))
    fn1 = f1[length(f1)]
    dir = paste(collapse="\\/", c(f1[1:(length(f1)-1)]) )

  USTA= unique(stns)
  nn =length(ma)
  pcol=rep(1, nn)
  if(length(USTA)==1)
    {
      pcol=rep(1, nn)
    }
  else
    {
      for(m in 1:length(USTA))
        {
          pcol[!is.na(match( stns, USTA[m]))] = 2+m
        }
    }
  ok = order(notes)
  wintim=  info$jd[1] + info$hr[1]/24+ info$mi[1]/(24*60)+(info$sec[1]+info$msec[1]/1000+info$t1[1]-info$off[1])/(24*3600)
  ex = seq(0,length(ascd[[1]])-1)*info$dt[1]

  dat=NULL
    
  ftime = paste(sep="_", info$yr,info$mo,info$dom,info$hr,info$mi,info$sec,info$msec)


  ###################  compnent names
  ###################         in order to be consistent with the pickfiles
  ################      for matching, here we convert all compnent names to a standard convention
  ###########   the convention is: V N E    NOT Z
  ############   usually 4,5,6 = VNE sometimes 1,2,3  = VNE
  ####   sometimes A = 1

 Orig.comps = comps
  comps = fixcomps(Orig.comps)
  
  GFIL = list(JSTR=ascd, STNS=stns, dir=dir, ifile=fn1, COMPS=comps, OCOMPS=Orig.comps,
    dt=dt, KNOTES=notes, info=info,dat=dat, nn=nn, ex=ex, pcol=pcol, ok=ok, wintim=wintim,  ftime=ftime, units=units )
  
  invisible(GFIL)
  
}

