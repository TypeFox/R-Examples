#' @export
#' @import lme4
arfimaMLM <-
function(formula, data, timevar
         , d="Hurst", arma=NULL
         , ecmformula=NULL, decm="Hurst"
         , drop=5, report.data=TRUE
         , ...){
  
  # extract variables from formula
  varlist <- extractVars(formula, ecmformula)
  
  # prepare dataset
  new <- arfimaPrep(data=data, timevar=timevar
                    , varlist.mean=varlist$mean, varlist.fd=varlist$fd
                    , varlist.xdif=varlist$xdif, varlist.ydif=varlist$ydif
                    , d=d, arma=arma, ecmformula=ecmformula, decm=decm, drop=drop)
  
  # estimate multilevel model  
  res <- lmer(formula, data=new$data.merged, ...)

  # output
  if (report.data==T){
    out <- list(result=res,d=new$d,arma=new$arma,ecm=new$ecm
                ,data.mean=new$data.mean,data.fd=new$data.fd,data.merged=new$data.merged)
  } else out <- list(result=res,d=new$d,arma=new$arma,ecm=new$ecm)
  class(out) <- "arfimaMLM"
  out
}
