fitQmap <- function(obs,mod,
                    method=c("PTF","DIST","RQUANT","QUANT","SSPLIN"),...){
  method <- match.arg(method)
  ffun <- match.fun(paste("fitQmap",method,sep=""))
  ffun(obs,mod,...)
}
