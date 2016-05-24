derv <-
function(func,x0,n=1,complex=FALSE,tol=.Machine$double.eps^0.5){
final<-NULL
for(i in 1:length(x0)){
if(!(is.numeric(x0[i])==TRUE | is.complex(x0[i])==TRUE)){stop("Evaluating the derivative at a non-numeric point")}
cauchy<-paste("function(",names(formals(func)),"){(",body2string(func),")/(",names(formals(func)),"-",x0[i],")^(",n,"+1)}",collapse="",sep="")
result<-disc.integrate(eval(parse(text=cauchy)),x0[i],rel.tol=tol)$value*factorial(n)/(2*pi*1i)
if(complex==TRUE){final<-c(final,result)}else{final<-c(final,Re(result))}
}
return(final)
}
