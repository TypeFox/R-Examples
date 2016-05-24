dc <- function(x,d,ch="&") { # d=0: "x&0"
  frac <- function(x,d) {  # fractional part
    res <- abs(x-trunc(x))
    if (!missing(d)) res <- round(10^d*res)
    res
  }
  d  <- max(d,1)
  fr <- frac(x,d)
#  fr <- ifelse (fr==0, paste(rep("0",d),collapse=""), fr)  to make it behave like dcn
  paste(trunc(x),ch,fr,sep="")
}

dcn <- function(x,d,ch="&") { # d=0: no "&" !!
  d  <- max(d,1)
  s <- sapply(x,function(x) eval(parse(text = paste("sprintf('%.", d, "f',", x, ")", sep = ""))))
  gsub( ".",ch, s )
}

mpf <- function(r,after) {paste(if (r<0) "-" else "+", eval(parse(text = paste("sprintf('%.",after,"f',abs(",r,"))",sep=""))))}
