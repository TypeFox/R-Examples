"estDesign" <-
function( n, smax, p.tr, biasrest=0.05)
{

if(length(n)!=1 || (n<=1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer greater than 1")}
if(length(p.tr)!=1 || p.tr>1 || p.tr<0){stop("true proportion p.tr must be specified as a single number between 0 and 1")}
if(length(smax)!=1 || (smax<=1 | abs(round(smax)-smax) > 1e-07)){stop("the maximal group size allowed in calculations must be a single integer greater than 1")}
if(length(biasrest)!=1 || biasrest>=1 || biasrest<0){stop("the maximally allowed bias(p) specified in biasrest must be a single number between 0 and 1, usually should be close to 0")}

 for (i in 2:smax)
  {
  temp<-msep(n=n, p.tr=p.tr,s=i)

   if(temp$bias > biasrest)
     {cat("maximal group size within bias restriction is s =",i-1,"\n")
       return(msep(n=n, p.tr=p.tr,s=i-1))} 

   if(i>=2 && temp$mse > msep(n=n, p.tr=p.tr,s=i-1)$mse)
     {cat("group size s with minimal mse(p) =",i-1,"\n")
      return(msep(n=n, p.tr=p.tr,s=i-1))}

   if (i==smax && temp$mse <= msep(n=n, p.tr=p.tr,s=i-1)$mse)
     {cat(" minimal mse(p) is achieved with group size s >= smax","\n")
      return(msep(n=n, p.tr=p.tr,s=i-1))}
  }
}

