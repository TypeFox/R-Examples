print.CAshiny<-
  function(x,...){
    res.shinyca=x
    if(!inherits(res.shinyca,"CAshiny"))
      stop("non convenient data")
    cat("Results for the CA with Factoshiny\n")
    cat("You can use it to fine your app the way you left it\n")
    cat("\n")
    cat("Corresponding script : \n")
    cat(res.shinyca$code1,"\n")
    cat("\n")
    cat(res.shinyca$code2,"\n")
  }