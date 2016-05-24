"Fisher1sided" <-
function(ctable,side="great"){

if(length(ctable)!=4) stop("2x2 table must have 4 entries")

### yobs and margins
yobs <- ctable[1]
n1 <- ctable[1]+ctable[2]
n2 <- ctable[1]+ctable[3]
nn <- sum(ctable)

########### probs and counts under the null

mx <- max(0,n1+n2-nn)
mn <- min(n1,n2)
ynull <- c(mx:mn)
dp <- dhyper(ynull,n2,nn-n2,n1) 

############ get p-value

dp <- round(dp,d=9)

ind.obs <- match(yobs,ynull)
dp.obs <- dp[ind.obs]

print(ynull)

if(side=="great") pval.Fish <- sum(dp[ind.obs:length(ynull)])
else if(side=="less") pval.Fish <- sum(dp[1:ind.obs])
else stop("side must be 'great' or 'less'")


return(list(pval.Fish=pval.Fish,Prob=dp.obs,yobs=yobs,n1=n1,n2=n2,nn=nn))

}

