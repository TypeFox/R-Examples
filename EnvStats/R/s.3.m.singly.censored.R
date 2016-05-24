s.3.m.singly.censored <-
function (r, alf, bet) 
{
    (r^2 * g1.m.singly.censored(r, 2.5, alf, bet) - g1.m.singly.censored(r, 
        1.5, alf, bet))/g1.m.singly.censored(r, 0.5, alf, bet)
}
