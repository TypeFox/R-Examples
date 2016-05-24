LASSO <-
function(x, y, family=c("gaussian","binomial","poisson","multinomial","cox"), 
			alpha = 1, nlambda = 100, lambda.min.ratio, standardize = TRUE, 
			size.max, a = 2, npermutation = 0, npermutation.max, min.nonsignificant.counts)
{
#
#a<- c(2, 1+log(length(y))/2, log(length(y)))
#IC: AIC, GIC, BIC

family<- match.arg(family, c("gaussian","binomial","poisson","multinomial","cox"))
alpha<- as.double(alpha)
nlambda<- as.integer(nlambda)
y<- drop(y)
np<- dim(x)
nobs<- as.integer(np[1])
nvars<- as.integer(np[2])
IC<- paste("GIC", round(a, digits=2), sep="")
IC[a == 2]<- "AIC"
IC[round(a, digits=2) == round(log(nobs), digits=2)]<- "BIC"


if (is.null(colnames(x))) colnames(x)<- paste("V", seq(nvars), sep = "")
if (missing(size.max)) size.max<- min(c(nobs-1, nvars))
if (missing(lambda.min.ratio)) lambda.min.ratio<- ifelse(nobs<nvars, 0.01, 0.0001)

if (missing(min.nonsignificant.counts)) min.nonsignificant.counts<- 10  #?10, ?20
if (missing(npermutation.max)) npermutation.max<- npermutation

#######
#(A) pvalue (nominal pvalue)
yp<- y
#(1) generate a set of models.
	fit<- glmnet(x, yp, family = family, alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, standardize = standardize)
	vselected<- (as.matrix(fit$beta)!=0)
	vselected<- unique(vselected, MARGIN=2) 
	vselected<- vselected[,colSums(vselected)<=size.max, drop = FALSE]
	colnames(vselected)<- colSums(vselected)

#(2) logLik and p-value
	#logL<- rep(NA, ncol(vselected))  #logL[i]<- logLik(fit)[1]
	pval<- rep(NA, ncol(vselected))
	#Fval<- rep(NA, ncol(vselected))
	ICval<- matrix(NA, nrow = ncol(vselected), ncol = length(a))
	colnames(ICval)<- IC
	for (i in 1:ncol(vselected)){
		xs<- x[,vselected[,i], drop=FALSE]
		if (family == "gaussian") {
			if (all(!vselected[,i])) {
			fit<- lm(yp~1)
			pval[i]<- 1
			} else {
			fit<- lm(yp~xs)
			pval[i]<- anova(fit)[1, "Pr(>F)"]
			}
			#ICval[i,]<- sapply(a, function(ak) extractAIC(fit, k = ak)[2])
		} else {
			if (all(!vselected[,i])) {
			fit<- glm(yp~1, family = family)
			pval[i]<- 1
			#ICval[i,]<- sapply(a, function(ak) extractAIC(fit, k = ak)[2])
			} else {
			fit<- glm(yp~xs, family = family)
			pval[i]<- 1 - pchisq(fit$null.deviance - fit$deviance, df = fit$df.null - fit$df.residual)  #LR testing
			#fit<- glmnet(xs, yp, family = family, lambda = 1e-5)
			#nulldev<- fit$nulldev
			#dev<- (1 - fit$dev.ratio)*nulldev
			#pval[i]<- 1 - pchisq(nulldev - dev, fit$df)
			#ICval[i,]<- -2*(dev - nulldev) + fit$df*a???
			}
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
	fit<- glmnet(x, yp, family = family, alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, standardize = standardize)
	vselected<- (as.matrix(fit$beta)!=0)
	vselected<- unique(vselected, MARGIN=2) 
	vselected<- vselected[,colSums(vselected)<=size.max, drop = FALSE]
	colnames(vselected)<- colSums(vselected)

#(2) logLik and p-value
	#logL<- rep(NA, ncol(vselected))  #logL[i]<- logLik(fit)[1]
	pval<- rep(NA, ncol(vselected))
	#Fval<- rep(NA, ncol(vselected))
	ICval<- matrix(NA, nrow = ncol(vselected), ncol = length(a))
	colnames(ICval)<- IC
	for (i in 1:ncol(vselected)){
		xs<- x[,vselected[,i], drop=FALSE]
		if (family == "gaussian") {
			if (all(!vselected[,i])) {
			fit<- lm(yp~1)
			pval[i]<- 1
			} else {
			fit<- lm(yp~xs)
			pval[i]<- anova(fit)[1, "Pr(>F)"]
			}
		} else {
			if (all(!vselected[,i])) {
			fit<- glm(yp~1, family = family)
			pval[i]<- 1
			} else {
			fit<- glm(yp~xs, family = family)
			pval[i]<- 1 - pchisq(fit$null.deviance - fit$deviance, df = fit$df.null - fit$df.residual)  #LR testing
			}
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

list(nonsignificant.counts = counts, pvalue.empirical = permpvalue, pvalue.nominal = NA, vs = testvs, total.permutation = jth, family = family)
}

