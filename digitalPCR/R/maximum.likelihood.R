maximum.likelihood <-
function(param, pos, neg, dilution) #p,
{
min.v=Inf
for (i in 1:nrow(param))
{
la =param[i,1]
tt =param[i,2]
diff=TL(pos, neg, la, tt, dilution)#total likelihood
if (is.na(diff)){
next
}
if (diff < min.v){
min.v=diff
sol=param[i,]
}
}
return (list(sol))
}
