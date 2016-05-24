Tstat <-
function(s1,s2,phi,n1,n2)
{
num = s2 - s1
mu = (s1+s2)/(n1+n2)
den = sqrt((n1+n2)*mu*(1+phi*mu))
return(abs(num/den))
}
