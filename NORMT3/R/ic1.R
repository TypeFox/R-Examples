"ic1" <-
function (p,d) 
{
cc <- complex(real=0, imaginary=p/2)

(sqrt(pi)/2)*exp(-(p^2)/4)*(1-Re(erf(cc) - erf(cc-d)))

}
