'binomest'<-function(x, ...){UseMethod("binomest")}

'binomest.default' <-
function(x, n, names=NULL, method="Wald", success=NULL, ...)
{

if(is.null(success))
 {success<-"success"}
else
 {success<-paste(as.character(success), collapse="")}

method<-match.arg(arg = method,
      choices = c("Wald", "ADD1", "ADD2"))

Ntot <- sum(n)
k <- length(n)


if(is.null(names))
 {names<-paste(1:k)}


if(length(x) != k || length(names)!= k)
 {stop("arguments x, n, and names must be vectors of equal length")}

Y<-x

estimate <- Y/n

if(method=="Wald")
{ estp <- estimate
  varp <- estp*(1-estp)/n

 which0<-which(Y==0)
 whichn<-which(Y==n)

 estcor <- estp
 estcor[which0] <- 0.5/n[which0]
 estcor[whichn] <- (n[whichn]-0.5)/n[whichn]

 varcor <- estcor*(1-estcor)/n
}

if(method=="ADD1")
{
  estp <- (Y+0.5)/(n+1)
  varp <- estp*(1-estp)/(n+1)
  varcor <- varp
}

if(method=="ADD2")
{
  estp <- (Y+1)/(n+2)
  varp <- estp*(1-estp)/(n+2)
  varcor <- varp
}


names(Y) <- names(n) <- names(estimate) <- names(estp) <- names(varp) <- names(varcor) <- names

return(list(
Y=Y,
n=n,
estimate=estimate,
estp=estp,
varp=varp,
varcor=varcor,
names=names,
success=success
))

}




'binomest.formula' <- function(formula, data, method="Wald", success=NULL, ...)
{

if(length(formula)!=3)
 {stop("formula must be of shape 'response ~ treatment'") }

mf<-model.frame(formula=formula, data=data)

mf[,2]<-as.factor(mf[,2])

mf[,1]<-as.factor(mf[,1])

tab<-table(mf)

 if(length(dimnames(tab)[[1]])!=2)
  {stop("The specified response variable should have two levels")}


if(is.null(success))
 {

  if(length(dimnames(tab)[[1]])!=2)
   {warning("The response variable has more than 2 levels. The first will be taken as success")}

  n<-colSums(tab)
  x<-tab[1,]
  names<-dimnames(tab)[[2]]
  success<-dimnames(tab)[[1]][1]
 }
else
 {
  success<-as.character(success)
  if(length(success)!=1)
   {stop("Argument success should be a single character string,
 naming one of the levels of the response variable")}

  if(all(dimnames(tab)==success))
   {stop(paste("Level", success, "could not be found in the specified response variable"))}

 n<-colSums(tab)
 x<-tab[success,]
 names<-dimnames(tab)[[2]]
 }


out <- binomest.default(x=x, n=n, names=names,
 method=method, success=success)

return(out)
}



'binomest.table' <-
function(x, method="Wald", success=NULL, ...)
{

if(length(dim(x))!=2)
 {stop("x must be a table with two dimensions")}

if(!any(dim(x)==2))
 {stop("None dimension of x equals 2, this function can only handle 2xk tables")}

respdim<-which(dim(x)==2)

groupdim<-which(dim(x)>=2)

if(respdim==1)
 {
  if(is.null(success))
   {
   resp <- as.numeric(x[1,])
   n<-as.numeric(colSums(x))
   names(resp) <- names(n) <- dimnames(x)[[2]]
   success <- dimnames(x)[[1]][1]
   }
  else
   {
   
   success<-as.character(success)
   if(length(success)!=1)
     {stop("success must be a single character string")}

   if(all(dimnames(x)[[1]]!=success))
     {stop(paste("success=", success,"can not be found in", names(dimnames(x)), collapse="" ))}

   resp <- as.numeric(x[success,])
   n<-as.numeric(colSums(x))
   names(resp) <- names(n) <- dimnames(x)[[2]]
   }
 }
 else
  {
  if(is.null(success))
   {
   resp <- as.numeric(x[,1])
   n <- as.numeric(rowSums(x))
   names(resp) <- names(n) <- dimnames(x)[[1]]
   success <- dimnames(x)[[2]][1]
   }
  else
   {
   
   success<-as.character(success)
   if(length(success)!=1)
     {stop("success must be a single character string")}

   if(all(dimnames(x)[[2]]!=success))
     {stop(paste("success=", success,"can not be found in", names(dimnames(x)), collapse="" ))}

   resp <- as.numeric(x[,success])
   n <- as.numeric(rowSums(x))
   names(resp) <- names(n) <- dimnames(x)[[1]]
   } 

  }

out <- binomest.default(x=resp, n=n, names=names(resp), method=method, success=success)

return(out)

}





