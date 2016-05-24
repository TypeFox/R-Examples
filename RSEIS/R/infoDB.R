infoDB<-function(DB)
  {

    if(any(DB$yr>2100))
      {

        print("Problems with DB....need to check")
        ww = which(DB$yr>2100)
        print(ww)


      }

    usta = unique(DB$sta)
    ucomp = unique(DB$comp)
    times = range(c(DB$t1, DB$t2))

    cat("Unique Stations:\n")
    cat(usta, sep=" "); cat("\n")
    cat("Unique Components:\n")
    cat(ucomp, sep=" "); cat("\n")
    cat("Times(jd):\n")
    cat(times, sep=" "); cat("\n")


    RD = rangedatetime(DB)
   cat("Dates:\n")
D1 = showdatetime(RD[[1]])
D2 = showdatetime(RD[[2]])

   ##  print(D1)
   ##  print(D2)
   
    invisible(list(usta=usta, ucomp=ucomp, start=D1, end=D2, at1=times[1], at2=times[2] ))


  }
