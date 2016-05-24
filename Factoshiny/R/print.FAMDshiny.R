print.FAMDshiny<-
  function(x,...){
    res.shinyfamd=x
    if(!inherits(res.shinyfamd,"FAMDshiny"))
      stop("non convenient data")
    cat("Results for the FAMD with Factoshiny\n")
    cat("You can use it to fine your app the way you left it\n")
    cat("\n")
    cat("Corresponding script : \n")
    cat(res.shinyfamd$code1,"\n")
    cat("\n")
    cat(res.shinyfamd$code2,"\n")
    cat(res.shinyfamd$code3,"\n")
    cat(res.shinyfamd$code4,"\n")
  }