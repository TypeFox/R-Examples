"LD2sided.pvals" <-
function(ctable){

if(length(ctable)!=4) stop("2x2 table must have 4 entries")

###################################
##### check first if there will be problem with delta p-values
n1 <- ctable[1]+ctable[2]
n2 <- ctable[1]+ctable[3]
nn <- sum(ctable)
mx <- max(0,n1+n2-nn)
mn <- min(n1,n2)
ynull <- c(mx:mn)
if(sum((nn-n1-n2+ynull)==0)>0){
 ctable <- c(ctable[3],ctable[4],ctable[1],ctable[2])
 print("relabelling table rows so delta well-behaved")
}
n1 <- ctable[1]+ctable[2]
n2 <- ctable[1]+ctable[3]
nn <- sum(ctable)
mx <- max(0,n1+n2-nn)
mn <- min(n1,n2)
ynull <- c(mx:mn)
if(sum((nn-n1-n2+ynull)==0)>0){
 ctable <- c(ctable[2],ctable[1],ctable[4],ctable[3])
 print("relabelling table columns so delta well-behaved")
}
###################################

### yobs and margins
yobs <- ctable[1]
n1 <- ctable[1]+ctable[2]
n2 <- ctable[1]+ctable[3]
nn <- sum(ctable)

meanY <- n1*n2/nn
modeY <- floor((n1+1)*(n2+1)/nn)

########### probs and counts under the null

mx <- max(0,n1+n2-nn)
mn <- min(n1,n2)
ynull <- c(mx:mn)
dp <- dhyper(ynull,n2,nn-n2,n1) 

########### different LD measures (for all possible tables with same margins)

xx <- ynull/nn
aa <- n1/nn
bb <- n2/nn

rho <- xx*(1+xx-aa-bb)/((aa-xx)*(bb-xx))
dd <- xx-aa*bb
dmax <- as.numeric(dd>=0)*min(aa*(1-bb),bb*(1-aa)) +
	  as.numeric(dd<0)*min(aa*bb,(1-aa)*(1-bb))
dprime <- dd/dmax
rr <- dd/sqrt(aa*(1-aa)*bb*(1-bb))
yuleQ <- 1-2/(rho+1)

#diff1 <- dd/(aa*(1-aa))
#diff2 <- dd/(bb*(1-bb))
#delta1 <- dd/(aa*(1+xx-aa-bb))
delta2 <- dd/(bb*(1+xx-aa-bb))

LR <- log( (ynull/(nn*aa*bb))^ynull ) +
	log( ((n1-ynull)/(nn*aa*(1-bb)))^(n1-ynull) ) +
	log( ((n2-ynull)/(nn*bb*(1-aa)))^(n2-ynull) ) +
	log( ((nn-n1-n2+ynull)/(nn*(1-aa)*(1-bb)))^(nn-n1-n2+ynull) )
LR <- 2*LR

############ get usual 2-sided p-values for the 6 stats

dp <- round(dp,d=9)
LR <- round(LR,d=9)
rr <- round(rr,d=9)
dprime <- round(dprime,d=9)
delta2 <- round(delta2,d=9)
yuleQ <- round(yuleQ,d=9)

 ind.obs <- match(yobs,ynull)

 dp.obs <- dp[ind.obs]
 LR.obs <- LR[ind.obs]
 rr.obs <- rr[ind.obs]
 dprime.obs <- dprime[ind.obs]
 delta2.obs <- delta2[ind.obs]
 yuleQ.obs <- yuleQ[ind.obs]

 pval.Fish <- sum(dp[dp<=dp.obs])
 pval.LR <- sum(dp[LR>=LR.obs])
 pval.r <- sum(dp[abs(rr)>=abs(rr.obs)])
 pval.Dprime <- sum(dp[abs(dprime)>=abs(dprime.obs)])
 pval.delta <- sum(dp[abs(delta2)>=abs(delta2.obs)])
 pval.yuleQ <- sum(dp[abs(yuleQ)>=abs(yuleQ.obs)])

############# conditional 2-sided p-values

weightleft <- sum(dp[ynull<=meanY])
weightright <- sum(dp[ynull>=meanY])

if(yobs==meanY) pval.cond <- 1
else if(yobs<meanY) pval.cond <- sum(dp[1:ind.obs])/weightleft
else pval.cond <- sum(dp[ind.obs:length(ynull)])/weightright

return(list(pval.cond=pval.cond,pval.Fish=pval.Fish,pval.LR=pval.LR,pval.r=pval.r,
	pval.Dprime=pval.Dprime,pval.delta=pval.delta,pval.Q=pval.yuleQ,
	Prob=dp.obs,LR=LR.obs,r=rr.obs,
	Dprime=dprime.obs,delta=delta2.obs,Q=yuleQ.obs,
	yobs=yobs,n1=n1,n2=n2,nn=nn))

}

