PLS <-
function(x, y, scale = FALSE, ncomp, varpercent, 
	npermutation = 100, npermutation.max, min.nonsignificant.counts)
{
#varpercent = 0.80
#always centralizing x and y
x<- scale(x, scale = scale)
y<- y-mean(y)

if (missing(min.nonsignificant.counts)) min.nonsignificant.counts<- 10  #?20
if (missing(npermutation.max)) npermutation.max<- npermutation
if (missing(ncomp) & missing(varpercent)) ncomp<- min(dim(x))
if (!missing(ncomp)) {
	stopifnot(ncomp >= 1)
	ncomp<- floor(ncomp)
	names(ncomp)<- paste("PLS", ncomp, sep="")
	}
else ncomp<- NULL
if (!missing(varpercent)) stopifnot(varpercent >0 & varpercent <= 1)


#(1) test score and pvalue (nominal pvalue)
yp<- y
if (missing(varpercent)) {
	ncompvar<- NA
	fit<- plsr(yp~x, scale = FALSE, ncomp =  max(ncomp))
	ncomptest<- ncomp
	}
else {
	fit<- plsr(yp~x, scale = FALSE)
	indx<- cumsum(explvar(fit)) >= (100*varpercent)
	ncompvar<- min(which(indx))
	names(ncompvar)<- paste("PLS", ncompvar, ".v", round(varpercent, digits=2), sep="")
	ncomptest<- c(ncomp, ncompvar)
	}
sco<- sapply(ncomptest, function(k) cor(yp, x%*%fit$coefficients[, , k]))
testscore<- abs(sco)
testpvalue<- NULL
#testscore<- sco*sqrt( (length(yp)-2)/(1-sco^2) ) #t-distribution with df of n-2. not correct??? cor(y, Ay)?
#testpvalue<- 2*pt(abs(testscore), df = length(yp)-2, lower.tail = FALSE)

#(2) permutation pvalue (empirical pvalue)
#permutation
permpvalue<- NULL
counts<- NULL
jth<- 0

if (npermutation >= 1) {
counts<- rep(0, length(ncomptest))
while ((jth < npermutation) | ((jth < npermutation.max) & (min(counts) < min.nonsignificant.counts)))
{
jth<- jth + 1
	yp <- y[sample(1:length(y), replace=FALSE)]
	if (missing(varpercent)) {
		fit<- plsr(yp~x, scale = FALSE, ncomp =  max(ncomp))
		ncomptest<- ncomp
		} else {
		fit<- plsr(yp~x, scale = FALSE)
		indx<- cumsum(explvar(fit)) >= (100*varpercent)
		ncomptest<- c(ncomp, min(which(indx)))
		}
	sco<- sapply(ncomptest, function(k) cor(yp, x%*%fit$coefficients[, , k]))
	permscore<- abs(sco)
	#permscore<- sco*sqrt( (length(yp)-2)/(1-sco^2) ) #t-distribution with df of n-2.
	counts<- counts + (permscore >= testscore)
	}#while
	permpvalue<- (1+ counts)/(1 + jth)
}#if

list(score = testscore, nonsignificant.counts = counts, pvalue.empirical = permpvalue, pvalue.nominal = testpvalue, total.permutation = jth, ncomp.varp = ncompvar)
}

