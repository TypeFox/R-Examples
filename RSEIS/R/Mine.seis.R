`Mine.seis` <-
  function(at1, at2, DB, grepsta, grepcomp,  kind=1, Iendian=1, BIGLONG=FALSE, CHOP=TRUE, verbose = FALSE, chtoken=NULL, statoken=NULL, RAW=FALSE)
{
  ##  find all the files that overlap the times times
#########   at1 and at2 should be single times in julian days
#######      if they are vectors, the first element will be used

##### grepsta, grepcomp should match what is on the file name, not
###  what is in the header of the seismic data


  ####  in some situations the chanel name and the station name are not embedded in the
  ####   file headers - in that case use the token from the file name

  if(missing(chtoken)) { chtoken=NULL }
  if(missing(statoken)) { statoken=NULL }
  
  if(missing(kind)) { kind = 1 }
  if(missing(CHOP)) { CHOP=TRUE }
  if(missing(verbose)) { verbose=FALSE }
  if(missing(Iendian)) { Iendian=1 }
     if(missing(BIGLONG)) { BIGLONG=FALSE }
 if(missing(RAW)) { RAW=FALSE  }

  

  if(is.null(DB$t1))
    {
      DB = T12.pix(DB)
    }

  
  ##if(at1[1]>365)
  ##  {
  ##    at1 
  ##
  ##  }
  
  w1 = which( at1[1]>=DB$t1 & at1[1]<DB$t2 )
  w2 = which( at2[1]>=DB$t1 & at2[1]<DB$t2 )
  w3 = which( DB$t1>=at1[1] & DB$t1<at2[1] )
  w4 = which( DB$t2>=at1[1] & DB$t2<at2[1] )

  if(length(c(w1, w2, w3, w4) )<1)
    {
      print("No time Match in DataBase")
      return()
    } 

  wi = unique(c(w1, w2, w3, w4))

  if(length(wi)<1) {
    print("No time  Match in DataBase")

    return()} 
  

  fn1 = DB$fn[wi]

  sta.comp = paste(DB$sta[wi], DB$comp[wi], sep=".")

  nsta = length(grepsta)

  if(is.null(grepcomp)) grepcomp= "*"
  
  ncomp = length(grepcomp)
  
  gi = vector()
  
  for(i in 1:nsta)
    {
      for(j in 1:ncomp)
        {
          gi = c(gi, grep(paste(sep=".", grepsta[i], grepcomp[j]), sta.comp ))
        }
    }
  
  if(length(gi)<1 ) {

    print("No grep(sta-comp) Match in DataBase, check that file names match DB")

    return() }
  
  fn2 = fn1[gi]
  
   if(verbose) print(fn2)

  
  KG4 = GET.seis(fn2, kind = kind,  Iendian=Iendian, BIGLONG= BIGLONG, HEADONLY=FALSE ,PLOT = -1, RAW=RAW)


 ## here we must insure that the channel names are unique

  
 if(verbose)
   { print(length(KG4) )
for(k in 1:length(KG4)) {
 cat( paste(KG4[[k]]$sta, KG4[[k]]$comp ), sep="\n")

}
     
   }

  

   if(verbose)
    {
      print("LENGTHS Retrieved:")
    for(ib in 1:length(KG4))
      {
        print(paste(ib, length(KG4[[ib]]$amp)))


      }
  }
 ## print("Mine GLUE.GET.seis")
  
  RR = GLUE.GET.seis(KG4)
 ##  print("Mine prepseis")
  GH=prepSEIS(RR)
  
  ##  the window is seconds from the begining of the traces
  ##   at1 and at2 are in julian days

  if(verbose)
    {
      Gdf = data.frame(
        yr=GH$info$yr,
        jd=GH$info$jd,
        mo=GH$info$mo,
        dom=GH$info$dom,
        hr=GH$info$hr,
        mi=GH$info$mi,
        sec=GH$info$sec,
        msec=GH$info$msec,
        dt=GH$info$dt,
        t1=GH$info$t1,
        t2=GH$info$t2,
        off=GH$info$off,
        n1=GH$info$n1,
        n2=GH$info$n2,
        n3=GH$info$n3,
        n=GH$info$n)

      print(Gdf)
    }


  if(is.null(attr(DB, "origyr")))  { origyr = min(DB$yr) }
  else {origyr=attr(DB, "origyr") }
  eday = EPOCHday(GH$info$yr, jd=GH$info$jd, origyr = origyr )

 #######  print("Mine after EPOCHday")
  
  ss1 = secdif(eday$jday, GH$info$hr, GH$info$mi, GH$info$sec-GH$info$off,
    at1,0, 0, 0)
  
  ss2 = secdif(eday$jday, GH$info$hr, GH$info$mi, GH$info$sec-GH$info$off,
    at2,0, 0, 0)
  
  
  win = c(min(ss1), max(ss2))


   if(win[1]<0) { win[1] = 0 }

  
  if(verbose)
    {

      print("Mine.seis internal eday=")
      print(eday)


      print("Mine.seis ss1=")
      
      print(ss1)

        print("Mine.seis ss2=")
      print(ss2)

          print("Mine.seis win=")
      
      print(win)
      print("GH")
      for(i in 1:length(GH$JSTR)) { print(paste(i, length(GH$JSTR[[i]]))) }
    }

  if(CHOP) {  HH = CHOP.SEISN(GH, WIN=win) }
  else
    {
      HH =GH

    }


  if(verbose)
    {
      print("HH")
      print(win)
      print(diff(win))
      
      for(i in 1:length(HH$JSTR)) { print(paste(i, length(HH$JSTR[[i]]))) }
    }

  return(HH)
  
}

