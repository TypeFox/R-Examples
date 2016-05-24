size.prop_test.two_sample=
function(n = NULL,p1,p2,alpha,power = NULL,alt)
# exactly one of "n" or "power" must be NULL
{if (is.null(n))
   {n=ceiling(power.prop.test(p1=p1, p2=p2,
    sig.level=alpha, power=power,
    alternative=alt)$n)
    n=ceiling((n/4)*(1+sqrt(1+4/(n*abs(p1-p2))))^2)
    return(n)
   }
 if (is.null(power))
   {q1=1-p1
   q2=1-p2
   p=(p1+p2)/2
   q=1-p
   c=p1*q1+p2*q2
   a=2*p*q
   delta=abs(p1-p2)
   b=sqrt(n)*delta
   if (alt == "two.sided")
       {alpha=alpha/2}
   power=pnorm((sqrt(a)*qnorm(1-alpha)-delta*sqrt(n))/sqrt(c))
   return(power)
   }
}


