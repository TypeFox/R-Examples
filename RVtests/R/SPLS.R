SPLS <-
function(x, y, scale=TRUE, ncomp, eta.grid, size.max, a = 2,
			npermutation = 0, npermutation.max, min.nonsignificant.counts)
{
#a<- c(2, 1+log(length(y))/2, log(length(y)))
IC<- paste("GIC", round(a, digits=2), sep="")
IC[a == 2]<- "AIC"
IC[round(a, digits=2) == round(log(length(y)), digits=2)]<- "BIC"

if (missing(ncomp)) ncomp<- 1
if (missing(eta.grid)) eta.grid<- seq(0.5, 0.9, by=0.1)
if (missing(size.max)) size.max<- min(c(ncol(x), length(y)-1))

if (missing(min.nonsignificant.counts)) min.nonsignificant.counts<- 10  #?10, ?20
if (missing(npermutation.max)) npermutation.max<- npermutation

#######
#(A) pvalue (nominal pvalue)
yp<- y
#(1) generate a set of models.
vselected<- NULL
for (eta in eta.grid){
	f<- spls(x, yp, eta=eta, K=ncomp, scale.x=scale, scale.y=FALSE)
	for (k in 1:ncomp){
		#vselected<- cbind(vselected, sapply(1:ncomp, function(k) f$betamat[[k]]!=0))
		betak<- (f$betamat[[k]]!=0)
		if (!any(is.na(betak))) vselected<- cbind(vselected, betak)
		}
	}
	vselected<- unique(vselected, MARGIN=2) 
	vselected<- vselected[,colSums(vselected)<=size.max, drop = FALSE]
	colnames(vselected)<- colSums(vselected)
	rownames(vselected)<- colnames(x)

#(2) logLik and p-value
	#logL<- rep(NA, ncol(vselected))  #logL[i]<- logLik(fit)[1]
	pval<- rep(NA, ncol(vselected))
	#Fval<- rep(NA, ncol(vselected))
	ICval<- matrix(NA, nrow = ncol(vselected), ncol = length(a))
	colnames(ICval)<- IC
	for (i in 1:ncol(vselected)){
		xs<- x[,vselected[,i], drop=FALSE]
		if (all(!vselected[,i])) {
			fit<- lm(yp~1)
			pval[i]<- 1
		} else {
			fit<- lm(yp~xs)
			pval[i]<- anova(fit)[1, "Pr(>F)"]
			}
		ICval[i,]<- sapply(a, function(ak) extractAIC(fit, k = ak)[2])
		}

# (3) selected by ICs: AIC, GIC, BIC  
	ks<- colSums(vselected)
	kIC<- apply(ICval, 2, which.min)
	pvalue<- pval[kIC]
	names(pvalue)<- IC
	vs<- vselected[,kIC, drop=FALSE]
	colnames(vs)<- IC
	
	testpvalue<- pvalue
	testvs<- vs
#######
#(B) permutation pvalue (empirical pvalue)
#permutation
permpvalue<- NULL
counts<- NULL
jth<- 0

if (npermutation >= 1) {
counts<- rep(0, length(a))
while ((jth < npermutation) | ((jth < npermutation.max) & (min(counts) < min.nonsignificant.counts)))
{
jth<- jth + 1
	yp <- y[sample(1:length(y), replace=FALSE)]
#(1) generate a set of models.
vselected<- NULL
for (eta in eta.grid){
	f<- spls(x, yp, eta=eta, K=ncomp, scale.x=scale, scale.y=FALSE)
	for (k in 1:ncomp){
		#vselected<- cbind(vselected, sapply(1:ncomp, function(k) f$betamat[[k]]!=0))
		betak<- (f$betamat[[k]]!=0)
		if (!any(is.na(betak))) vselected<- cbind(vselected, betak)
		}
	}
	vselected<- unique(vselected, MARGIN=2) 
	vselected<- vselected[,colSums(vselected)<=size.max, drop = FALSE]
	colnames(vselected)<- colSums(vselected)
	rownames(vselected)<- colnames(x)

#(2) logLik and p-value
	#logL<- rep(NA, ncol(vselected))  #logL[i]<- logLik(fit)[1]
	pval<- rep(NA, ncol(vselected))
	#Fval<- rep(NA, ncol(vselected))
	ICval<- matrix(NA, nrow = ncol(vselected), ncol = length(a))
	colnames(ICval)<- IC
	for (i in 1:ncol(vselected)){
		xs<- x[,vselected[,i], drop=FALSE]
		if (all(!vselected[,i])) {
			fit<- lm(yp~1)
			pval[i]<- 1
		} else {
			fit<- lm(yp~xs)
			pval[i]<- anova(fit)[1, "Pr(>F)"]
			}
		ICval[i,]<- sapply(a, function(ak) extractAIC(fit, k = ak)[2])
		}

# (3) selected by ICs: AIC, GIC, BIC  
	ks<- colSums(vselected)
	kIC<- apply(ICval, 2, which.min)
	pvalue<- pval[kIC]
	names(pvalue)<- IC
	vs<- vselected[,kIC, drop=FALSE]
	colnames(vs)<- IC

	counts<- counts + (testpvalue >= pvalue)
	}#while
	permpvalue<- (1+ counts)/(1 + jth)
}#if

list(nonsignificant.counts = counts, pvalue.empirical = permpvalue, pvalue.nominal = NA, vs = testvs, total.permutation = jth)
}

