g2.m.singly.censored <-
function (r, kk, alf, bet) 
{
    gamma(alf + kk)/(2 * gamma(alf) * bet^kk) * (sign(r) * pbeta(r^2/(2 * 
        bet + r^2), 0.5, alf + kk) + 1)
}
