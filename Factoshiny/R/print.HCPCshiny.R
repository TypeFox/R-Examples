print.HCPCshiny<-
  function(x,...){
    res.shinyhcpc=x
    if(!inherits(res.shinyhcpc,"HCPCshiny"))
      stop("non convenient data")
    cat("Results for the HCPC with Factoshiny\n")
    cat("You can use it to fine your app the way you left it\n")
    cat("\n")
    cat("Corresponding script : \n")
    print(res.shinyhcpc$anafact)
    cat(res.shinyhcpc$code1,"\n")
    cat("\n")
    cat(res.shinyhcpc$code2,"\n")
    cat(res.shinyhcpc$code3,"\n")
    cat(res.shinyhcpc$code4,"\n")
  }