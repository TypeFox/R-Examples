"fun.fmkl.nr" <-
function(x, l1, l2, l3, l4, init, error = 1e-010)
{
u <- init
if(l3 == 0) {
if(l4 == 0) {
quants <- l1 + (log(u) - log(1 - u))/l2
}
else {
quants <- l1 + (log(u) - ((1 - u)^l4 - 1)/l4)/l2
}
}
else {
if(l4 == 0) {
quants <- l1 + ((u^l3 - 1)/l3 - log(1 - u))/l2
}
else {
quants <- l1 + ((u^l3 - 1)/l3 - ((1 - u)^l4 - 1)/l4)/l2
}
}
quants
}

