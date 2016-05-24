info.seis<-function(GH)
  {

    A = dateStamp(GH$info)
    B = paste(sep=" ", 1:length(GH$STNS), GH$STNS, GH$COMPS, GH$dt, GH$info$n )

    B2 = paste(sep=" ",B, A)
    

    cat(B2, sep="\n")

    

##########      info.seis(GH)
    
  }
