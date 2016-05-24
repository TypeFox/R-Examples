g1.m.singly.censored <-
function (r, kk, alf, bet) 
{
    (((1/sqrt(2 * pi)) * bet^alf)/(bet + r^2/2)^(alf + kk) * 
        gamma(alf + kk))/gamma(alf)
}
