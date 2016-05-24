likelihood <-
function(po, ne, lamda, thres)
{
pneg=stats::ppois(thres-1, lamda)
ppos=1-pneg
return (log(ppos)*po+log(pneg)*ne) 
}
