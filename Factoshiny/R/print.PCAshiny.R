print.PCAshiny<-
function(x,...){
  res.shinypca=x
  if(!inherits(res.shinypca,"PCAshiny"))
    stop("non convenient data")
  cat("Results for the PCA with Factoshiny\n")
  cat("You can use it to fine your app the way you left it\n")
  cat("\n")
  cat("Corresponding script : \n")
  cat(res.shinypca$code1,"\n")
  cat("\n")
  cat(res.shinypca$code2,"\n")
  cat(res.shinypca$code3,"\n")
}