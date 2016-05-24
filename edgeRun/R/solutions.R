solutions <-
function(T0,phi,n1,n2,s2,piece)
{
w1 = 1 + phi*{T0^2}/(n1+n2)
w2 = 1 - phi*{T0^2}/(n1+n2)
a = w2
b = -2*s2*w1 - {T0^2}
c = w2*{s2^2} - {T0^2}*s2
if(piece==1) {result.lower = pmax(ceiling(pmin({-b-sqrt({b^2}-4*a*c)}/{2*a},{-b+sqrt({b^2}-4*a*c)}/{2*a})),0)
return(result.lower)}
if(piece==2) {result.upper = pmax(floor(pmax({-b-sqrt({b^2}-4*a*c)}/{2*a},{-b+sqrt({b^2}-4*a*c)}/{2*a})),0)
return(result.upper)}
}
