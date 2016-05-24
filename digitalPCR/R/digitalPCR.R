digitalPCR <-
function(pos, neg, dilution, Nboot, single.copy = c("FALSE", "TRUE"), upper.copy=100)
{
single.copy <- match.arg(single.copy)
#single.copy=FALSE
#print(single.copy)
tab=poisson.param(single.copy, upper.copy)

boot.result=NULL
for (i in 1:Nboot){
bsample=bootstrapping(pos, neg)
best=maximum.likelihood(tab, bsample[[1]], bsample[[2]], dilution)[[1]]
boot.result=rbind(boot.result, best)
}
sdcopy=round(stats::sd(boot.result[,1]),2) #standard deviation copy number 
meancopy=round(mean(boot.result[,1]),2) #mean copy number
print (paste('copy number', meancopy, 'SD', sdcopy))
return (list('mean copy number'=meancopy, 
'sd copy number'=sdcopy,
'copy numbers'=boot.result[,1], 
'thresholds'=boot.result[,2]))
}
