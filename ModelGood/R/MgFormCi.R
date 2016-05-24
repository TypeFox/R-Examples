MgFormCi <- function(lower,upper,brace=c("(",")"),sep=";",degenerated="--",digits=3,percent=TRUE,space=""){
  lower <- round(lower,digits)
  upper <- round(upper,digits)
  if (length(brace)==1) brace <- c(brace,brace)
  if (percent==TRUE) {
    lower <- round(100*lower,digits-2)
    upper <- round(100*upper,digits-2)
  }
  else{
    lower <- round(lower,digits)
    upper <- round(upper,digits)
  }
  out <- paste(brace[1],lower,sep,upper,brace[2],sep=space)
  ##   out <- switch(as.character(style),"1"=c("lower"=lower,"upper"=upper),"2"=paste("[",lower,";",upper,"]",sep=sep),"3"=paste("(",lower,";",upper,")",sep=sep),"4"=paste(lower,"-",upper,sep=sep),"5"=paste("CI-95%: ", lower,"-",upper,sep=sep))
  ##   if (style>1 && is.character(degenerated)) out[lower==upper] <- degenerated
  out
}
