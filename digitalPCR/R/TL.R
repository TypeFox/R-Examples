TL <-
function(pos, neg, lamda, thres, dilution)
{
rval =0
i = 1
for (dil in dilution){
rval = rval+ likelihood(pos[i], neg[i], lamda/dil, thres)
i = i+1
}
return (-rval)
}
