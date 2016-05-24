print.MCAshiny<-
  function(x,...){
    res.shinymca=x
    if(!inherits(res.shinymca,"MCAshiny"))
      stop("non convenient data")
    cat("Results for the MCA with Factoshiny\n")
    cat("You can use it to fine your app the way you left it\n")
    cat("\n")
    cat("Corresponding script : \n")
    cat(res.shinymca$code1,"\n")
    cat("\n")
    cat(res.shinymca$code2,"\n")
    cat(res.shinymca$code3,"\n")
    if(!is.null(res.shinymca$code4)){
      cat(res.shinymca$code4,"\n")
    }
  }