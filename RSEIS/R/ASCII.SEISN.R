
ASCII.SEISN<-function(GH, sel=1,  HEAD=TRUE  )
{
  if(missing(sel)) { sel = 1:length(GH$STNS) }
  if(missing(HEAD)) { HEAD=TRUE }


  Aname =  Zdate(GH$info, sel, 0)


  for(i in 1:length(sel))
    {
      j = sel[i]
      fn = paste(sep=".", Aname[1], GH$STNS[j], GH$COMPS[j])

      sec = GH$info$sec[i]
      hed = paste(GH$info$yr[i], GH$info$jd[i], GH$info$hr[i], GH$info$mi[i],  sec, GH$info$dt[i])
      
      cat(paste(i, fn, sep=" "), sep="\n")
      ACON = file(description =fn, open = "w")

      if(HEAD) { cat(hed, file = ACON, sep = "\n") }
         cat(GH$JSTR[[i]], file = ACON, sep = "\n")


         
         close(ACON)
         
       }



    }
