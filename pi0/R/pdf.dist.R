pdf.dist=function(f1, f2, method=c('Hellinger','abdif'))
{
    method=match.arg(method)
    if(method=='abdif') {
        f=function(x)abs(f1(x)-f2(x))
        ans=try( integrate(f, -Inf, Inf)$value )
        return( c(abdif=if(is.numeric(ans)) ans else NA_real_ ))
    }
    if(method=='Hellinger') {
        f=function(x) (sqrt(f1(x))-sqrt(f2(x)))^2
        ans=try( integrate(f, -Inf, Inf)$value )
         return( c(Hellinger=if(is.numeric(ans)) sqrt(.5*ans) else NA_real_))
   }
}
