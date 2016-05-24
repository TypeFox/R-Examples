### R code from vignette source 'sdc_guidelines.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sdc_guidelines.rnw:191-228
###################################################
#require(sdcMicro)
#load("../Daten/ses.RData")	
f1 <- function(x){
	truth <- weighted.mean(x$earningsMonth, x$GrossingUpFactor.x)
	SEQ <- seq(10,100,10)
	risk <- risk2 <- utility <- utility2 <- perturbed <- perturbed2 <- numeric(length(SEQ))
	j <- 0
	for(i in SEQ){
		j=j+1
		ad <- addNoise(x[,c("earnings","earningsMonth")], noise=i, method="restr")
		ad2 <- microaggregation(x[,c("earnings","earningsMonth")], aggr=j+1, method="pca")
		perturbed[j] <- weighted.mean(ad$xm[,2], x$GrossingUpFactor.x)
		perturbed2[j] <- weighted.mean(ad2$mx[,2], x$GrossingUpFactor.x)
		utility[j] <- dUtility(ad$x, ad$xm)		
		risk[j] <- dRisk(ad$x, ad$xm, k=0.01)
		utility2[j] <- dUtility(ad$x, ad2$mx)		
		risk2[j] <- dRisk(ad$x, ad2$mx, k=0.01)
	}	
	list(truth=truth, perturbed=perturbed, utility=utility, risk=risk, 
			perturbed2=perturbed2, utility2=utility2, risk2=risk2, SEQ=SEQ)
}
#set.seed(123)
#res <- f1(x)
#save(res, file="res.RData")
load("res.RData")
par(cex.lab=1.5, mar=c(5,4.5,1,0.1))
plot(cbind(res$risk, res$utility), type="l", 
	xlab="disclosure risk", ylab="information loss",
	xlim=c(0.08,0.26), ylim=c(0.1,1.95))
lines(cbind(res$risk2, res$utility2), lty=2)
text(x=res$risk, y=res$utility, res$SEQ)
text(x=res$risk2, y=res$utility2, 2:11)
text(x=0.22,y=0.5, "disclosive", cex=1.5)
text(x=0.21,y=1.8, "disclosive and worst data", cex=1.5)
text(x=0.1,y=0.5, "good", cex=1.5)
text(x=0.11,y=1.8, "worst data", cex=1.5)
legend("right", legend=c("method1","method2"), lty=c(1,2))	


###################################################
### code chunk number 2: freq
###################################################
#require(devtools)
require(sdcMicro)
require(xtable)
data(francdat)   ## toy data set
sdc <- createSdcObj(francdat, keyVars=c('Key1','Key2','Key3','Key4'), numVars=c('Num1','Num2','Num3'), w='w')
df <- cbind(francdat[,c(2,4,5,6,8)], get.sdcMicroObj(sdc, "risk")$individual)	
df$Key3[df$Key3==5] <- 2

colnames(df)[1:4] <- c("Age", "Location", "Sex", "Education")

#colnames(df)[ncol(df)] <- expression(hat(F)[k])
df <- xtable(df, digits=c(0,0,0,0,0,1,3,0,1), align = "|l|llll|l|l|ll|",
		caption="Example of sample and estimated population frequency counts.", 
		label="listingFreq")


###################################################
### code chunk number 3: freqprint
###################################################
print(df,include.rownames = getOption("xtable.include.rownames", TRUE), caption.placement="top")


###################################################
### code chunk number 4: sdc_guidelines.rnw:718-723
###################################################
data(testdata)
sdc <- createSdcObj(testdata,
		keyVars=c('urbrur','roof','walls','water','electcon','relat','sex'),
		numVars=c('expend','income','savings'), w='sampling_weight', hhId ='ori_hid')
print(sdc, "risk")


###################################################
### code chunk number 5: microaggregation
###################################################
df <- francdat[,c(1,3,7)]	
df <- cbind(df, microaggregation(df, aggr=2)$mx)
colnames(df)[4:6] <- paste("Mic",1:3, sep="")
df <- xtable(df, digits=c(0,2,3,0,2,2,1), align = "|l|lll|lll|",
	caption="Example of micro-aggregation. Columns 1-3 contain the original variables, columns 4-6 the micro-aggregated values.", 
	label="listingMicroaggregation")


###################################################
### code chunk number 6: allprint
###################################################
print(df,include.rownames = getOption("xtable.include.rownames", TRUE), caption.placement="top")


###################################################
### code chunk number 7: sdc_guidelines.rnw:1385-1387
###################################################
require(laeken, quiet=TRUE)
data(ses)


