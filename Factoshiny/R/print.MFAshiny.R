print.MFAshiny<-
  function(x,...){
    res.shinymfa=x
    if(!inherits(res.shinymfa,"MFAshiny"))
      stop("non convenient data")
    cat("Results for the MFA with Factoshiny\n")
    cat("You can use it to fine your app the way you left it\n")
    cat("\n")
    cat("Corresponding script : \n")
    print(res.shinymfa$ligne)
    cat("\n")
    cat(res.shinymfa$code1,"\n")
    cat(res.shinymfa$code2,"\n")
    if(!is.null(res.shinymfa$code3)){
      cat(res.shinymfa$code3,"\n")
    }
    cat(res.shinymfa$code4,"\n")
    if(!is.null(res.shinymfa$code5)){
      cat(res.shinymfa$code5,"\n")
    }
  }