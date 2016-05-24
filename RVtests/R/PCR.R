PCR <-
function(x, y, scale = FALSE, ncomp, varpercent, 
	npermutation = 100, npermutation.max, min.nonsignificant.counts)
{
#x=UDV': svd(x)
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
	names(ncomp)<- paste("PCR", ncomp, sep="")
} else ncomp<- NULL

UDV<- svd(x)
	u<- UDV$u
	d<- UDV$d

ncompvar<- NA
if (!missing(varpercent)) {
	stopifnot(varpercent >0 & varpercent <= 1)
	ncompvar<- min(which( cumsum(d^2) /sum(d^2) >= varpercent)) 
	names(ncompvar)<- paste("PCR", ncompvar, ".v", round(varpercent, digits=2), sep="")
	ncomp<- c(ncomp, ncompvar)
	}

u<- u[,1:max(ncomp),drop=FALSE]

#(1) test score and pvalue (nominal pvalue)
yp<- y
uy<- t(u)%*%yp
sco<- sapply(ncomp, function(k) cor(yp, u[,1:k,drop=FALSE]%*%uy[1:k]))
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
counts<- rep(0, length(ncomp))
while ((jth < npermutation) | ((jth < npermutation.max) & (min(counts) < min.nonsignificant.counts)))
{
jth<- jth + 1
	yp <- y[sample(1:length(y), replace=FALSE)]
	uy<- t(u)%*%yp
	sco<- sapply(ncomp, function(k) cor(yp, u[,1:k,drop=FALSE]%*%uy[1:k]))
	permscore<- abs(sco)
	#permscore<- sco*sqrt( (length(yp)-2)/(1-sco^2) ) #t-distribution with df of n-2.
	counts<- counts + (permscore >= testscore)
	}#while
	permpvalue<- (1+ counts)/(1 + jth)
}#if

list(score = testscore, nonsignificant.counts = counts, pvalue.empirical = permpvalue, pvalue.nominal = testpvalue, total.permutation = jth, ncomp.varp = ncompvar)
}

