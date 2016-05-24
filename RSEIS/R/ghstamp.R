ghstamp <-function(GH, sel, WIN=c(485,600) )
  {
    ftime  <-  Zdate(GH$info, sel, t1=0)
    p1 = paste(ftime , GH$STNS[sel], GH$COMPS[sel], WIN[1], WIN[2])
    return(p1)
  }
