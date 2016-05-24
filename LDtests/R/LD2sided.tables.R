"LD2sided.tables" <-
function(ctable){

if(length(ctable)!=4) stop("2x2 table must have 4 entries")

###yobs <- ctable[1]
n1 <- ctable[1]+ctable[2]
n2 <- ctable[1]+ctable[3]
nn <- sum(ctable)

meanY <- n1*n2/nn
#modeY <- floor((n1+1)*(n2+1)/nn)

########### probs and counts under the null

mx <- max(0,n1+n2-nn)
mn <- min(n1,n2)
ynull <- c(mx:mn)
dp <- dhyper(ynull,n2,nn-n2,n1) 

weightleft <- sum(dp[ynull<=meanY])
weightright <- sum(dp[ynull>=meanY])

xx <- ynull/nn
aa <- n1/nn
bb <- n2/nn

y12 <- n1-ynull
y21 <- n2-ynull
y22 <- nn-n1-n2+ynull

########### different LD measures (for all possible values of y)

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

xtable <- matrix(rep(0,length(ynull)*10),ncol=10)

xtable[,1] <- ynull
xtable[,2] <- y12
xtable[,3] <- y21
xtable[,4] <- y22
xtable[,5] <- dp
xtable[,6] <- LR
xtable[,7] <- rr
xtable[,8] <- dprime
xtable[,9] <- delta2
xtable[,10] <- yuleQ

xtable <- data.frame(xtable)
names(xtable) <- c("n11","n12","n21","n22","Prob","LR","r","D'","delta","Q")

############ get orderings 

ord1 <- order(-round(dp,d=9))
ord2 <- order(round(LR,d=9))
ord3 <- order(abs(round(rr,d=9)))
ord4 <- order(abs(round(dprime,d=9)))
ord5 <- order(abs(round(delta2,d=9)))
ord6 <- order(abs(round(yuleQ,d=9)))

rank1 <- rank(-round(dp,d=9))
rank2 <- rank(round(LR,d=9))
rank3 <- rank(abs(round(rr,d=9)))
rank4 <- rank(abs(round(dprime,d=9)))
rank5 <- rank(abs(round(delta2,d=9)))
rank6 <- rank(abs(round(yuleQ,d=9)))

##### practice tables
#ytable <- matrix(rep(0,length(ynull)*6),nrow=6)
#ytable[1,] <- ynull[ord1]
#ytable[2,] <- ynull[ord2]
#ytable[3,] <- ynull[ord3]
#ytable[4,] <- ynull[ord4]
#ytable[5,] <- ynull[ord5]
#ytable[6,] <- ynull[ord6]
#row.names(ytable) <- c("Fish","LR","|r|","|D'|","|delta|","|Q|")
#print(ytable)
ytable <- matrix(rep(0,length(ynull)*6),nrow=6)
ytable[1,] <- rank1[ord1]
ytable[2,] <- rank2[ord2]
ytable[3,] <- rank3[ord3]
ytable[4,] <- rank4[ord4]
ytable[5,] <- rank5[ord5]
ytable[6,] <- rank6[ord6]
row.names(ytable) <- c("Fish","LR","|r|","|D'|","|delta|","|Q|")
#print(ytable)

ytable <- matrix(rep(0,length(ynull)*6),nrow=6)

dum11 <- gsub("FALSE"," ",gsub("TRUE",")",duplicated(rank1[ord1])))
dum21 <- gsub("FALSE"," ",gsub("TRUE",")",duplicated(rank2[ord2])))
dum31 <- gsub("FALSE"," ",gsub("TRUE",")",duplicated(rank3[ord3])))
dum41 <- gsub("FALSE"," ",gsub("TRUE",")",duplicated(rank4[ord4])))
dum51 <- gsub("FALSE"," ",gsub("TRUE",")",duplicated(rank5[ord5])))
dum61 <- gsub("FALSE"," ",gsub("TRUE",")",duplicated(rank6[ord6])))

dum12 <- rep("",length(ynull))
dum22 <- rep("",length(ynull))
dum32 <- rep("",length(ynull))
dum42 <- rep("",length(ynull))
dum52 <- rep("",length(ynull))
dum62 <- rep("",length(ynull))

dum12[(1:length(ynull))[dum11==")"]-1] <- "("
dum22[(1:length(ynull))[dum21==")"]-1] <- "("
dum32[(1:length(ynull))[dum31==")"]-1] <- "("
dum42[(1:length(ynull))[dum41==")"]-1] <- "("
dum52[(1:length(ynull))[dum51==")"]-1] <- "("
dum62[(1:length(ynull))[dum61==")"]-1] <- "("

ytable[1,] <- paste(dum12,ynull[ord1],dum11,sep="")
ytable[2,] <- paste(dum22,ynull[ord2],dum21,sep="")
ytable[3,] <- paste(dum32,ynull[ord3],dum31,sep="")
ytable[4,] <- paste(dum42,ynull[ord4],dum41,sep="")
ytable[5,] <- paste(dum52,ynull[ord5],dum51,sep="")
ytable[6,] <- paste(dum62,ynull[ord6],dum61,sep="")

ytable <- data.frame(ytable)
row.names(ytable) <- c("Fish","LR","|r|","|D'|","|delta|","|Q|")

####### cat fn not necessary
####cat(t(ytable),fill=24)

if( max(table(rank1[ord1]))>2 ) 
	print("Brackets not displayed correctly for ties of more than 3")
if( max(table(rank2[ord2]))>2 ) 
	print("Brackets not displayed correctly for ties of more than 3")
if( max(table(rank3[ord3]))>2 ) 
	print("Brackets not displayed correctly for ties of more than 3")
if( max(table(rank4[ord4]))>2 ) 
	print("Brackets not displayed correctly for ties of more than 3")
if( max(table(rank5[ord5]))>2 ) 
	print("Brackets not displayed correctly for ties of more than 3")
if( max(table(rank6[ord6]))>2 ) 
	print("Brackets not displayed correctly for ties of more than 3")

############ get p-values

dp <- round(dp,d=9)
LR <- round(LR,d=9)
rr <- round(rr,d=9)
dprime <- round(dprime,d=9)
delta2 <- round(delta2,d=9)
yuleQ <- round(yuleQ,d=9)

ztable <- matrix(rep(0,length(ynull)*11),ncol=11)

for(i in 1:length(ynull)){

 yobs <- ynull[i]
 y12obs <- y12[i]
 y21obs <- y21[i]
 y22obs <- y22[i]

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

 if(yobs==meanY) pval.cond <- 1
 else if(yobs<meanY) pval.cond <- sum(dp[1:ind.obs])/weightleft
 else pval.cond <- sum(dp[ind.obs:length(ynull)])/weightright

 ztable[i,1] <- yobs
 ztable[i,2] <- y12obs
 ztable[i,3] <- y21obs
 ztable[i,4] <- y22obs
 ztable[i,5] <- pval.Fish
 ztable[i,6] <- pval.LR
 ztable[i,7] <- pval.r
 ztable[i,8] <- pval.Dprime
 ztable[i,9] <- pval.delta
 ztable[i,10] <- pval.yuleQ
 ztable[i,11] <- pval.cond 

}

ztable <- data.frame(ztable)
names(ztable) <- c("n11 obs","n12 obs","n21 obs","n22 obs","pval.Fish","pval.LR",
     "pval.r","pval.Dprime","pval.delta","pval.Q","pval.cond")

#################### print results

cat("\n LD measures for all possible 2x2 tables:",fill=T)
print(round(xtable,d=3))

cat("\n Orderings of 2x2 tables for different LD measures:",fill=T)
print(ytable)

cat("\n P-values for all possible 2x2 tables:",fill=T)
print(round(ztable,d=3))

}

