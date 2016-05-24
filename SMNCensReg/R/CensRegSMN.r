CensReg.SMN <- function(cc, x, y, LS=NULL, nu=3, delta=NULL, cens="left", dist="T", show.envelope="FALSE", error=0.0001,iter.max=300)
{
  if(cens=="left"){cens = 1}
  if(cens=="right"){cens = 2}
  if(cens=="interval"){cens = 3}
  type = dist 
  
  namesx <- ('x1     ')
  if(ncol(as.matrix(x))>1)
  {
  	for(i in 2:ncol(as.matrix(x))){namesx <- cbind(namesx, paste("x",i,"     ",sep=""))}
  }
	if(ncol(as.matrix(y)) > 1) stop("Only univariate linear regression supported!")
  if(ncol(as.matrix(cc)) > 1) stop("Only univariate linear regression supported!")
	if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
	if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
	if( (type != "T") && (type != "Normal") && (type != "PearsonVII") && (type != "Slash") && (type != "NormalC")) stop("Distribution family not supported. Check documentation!")

	for(i in 1:length(y))
  {
	  if(cens=="3")
    {
		  if(is.na(y[i])==FALSE && is.na(LS[i])==FALSE)
      {
			  if( (y[i]==-Inf) | (y[i]==Inf) | (LS[i]==Inf) | (LS[i]==-Inf) ) stop("This package does not support mixed types of censoring. For left censoring use cens=1. For right censoring use cens=2. For interval censoring use cens=3.")
			}
		}else
    {
			  if( (y[i]==-Inf) | (y[i]==Inf) ) stop("This package does not support mixed types of censoring. For left censoring use cens=1. For right censoring use cens=2. For interval censoring use cens=3.")
		}
  }	
	
	if( type == "T" | type == "Slash" )
  {
    if(length(nu) > 1) stop("nu parameter must be a scalar")
  	if(length(nu) == 0) stop("nu parameter must be provided.")
  	if(nu <= 0) stop("nu parameter must be positive.")
  }
  if( type == "PearsonVII" )
  {
  	if(length(nu) > 1) stop("nu parameter must be a scalar")
  	if(length(nu) == 0) stop("initial value for nu parameter must be provided in case of Pearson VII distribution.")
  	if(nu <= 0) stop("nu parameter must be positive.")
  	if(length(delta) > 1) stop("delta parameter must be a scalar")
  	if(length(delta) == 0) stop("delta parameter must be provided in case of Pearson VII distribution.")
  	if(delta <= 0) stop("delta parameter must be positive.")
  }
  if(type == "NormalC")
  {
    if(length(nu) !=2) stop("nu must be a bidimensional vector in case of Contaminated Normal distribution")
    if(nu[1] <=0 || nu[1] >= 1) stop("nu[1] must lies in (0,1)")
    if(nu[2] <=0 || nu[2] >= 1) stop("nu[2] must lies in (0,1)")
  }
  if( (cens != "1") && (cens != "2") && (cens != "3"))  stop("Censored type not supported. 1 for left censoring, 2 for right censoring and 3 for intervalar censoring.")
  
  out <- EM.Cens.Int(cc,x,y,LS,nu,delta,cens,type,error,iter.max)
  SE <- t(out$SE)
  SE <- round(t(SE),digits=5)
  param <- round(cbind(rbind(out$betas,out$sigma2),SE),digits=5)
  namespar <- colnames(x) 
  colx <- ncol(as.matrix(x))
  if(length(namespar)==0)namespar <- namesx[1:colx]
  dimnames(param) <- list(c(namespar,expression(sigma^2)),c("Estimates", "SE"))
  if( type=="PearsonVII")
  {
    sig2 <- round(t(cbind(out$nu,out$delta )),digits=5)
    dimnames(sig2) <- list(c(expression(nu),expression(delta)),"")
  }
  if( (type=="T") || (type=="Slash"))
  {
    nu1 <- matrix(round(out$nu,digits=5),ncol=1,nrow=1)
    row.names(nu1) <- "nu"
    colnames(nu1) <- " "
  }
  if( type=="NormalC")
  {
    nuf <- t(as.matrix(out$nu))
    nuf <- t(nuf)
    sig2 <- round(nuf,digits=5)
    dimnames(sig2) <- list(c(expression(nu1),expression(nu2)),"")
  }
  cat('\n') 
  cat('-------------------------------------------\n')
  cat('         EM estimates and SE \n')
  cat('-------------------------------------------\n')
  print(param)
  if(type!="Normal")
  {
    if(type=="T"|type=="Slash")
    {
      print(nu1)
    }
    else
    {
      print(sig2)
    }
  }
  cat('------------------------------------------\n')
  cat('\r \n')
  critFin <- c(out$logver, out$AIC, out$BIC, out$EDC)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","EDC"))
  cat('\n') 
  cat('Model selection criteria\n')
  cat('-------------------------------------------\n')
  print(critFin)
  cat('-------------------------------------------\n')
  cat('\r \n')
  if(show.envelope=="TRUE")
  {
    envelop <- EnvelopeRMT(cc,x,y,LS,nu,delta,cens=cens,type=type)
  }
  out
}
