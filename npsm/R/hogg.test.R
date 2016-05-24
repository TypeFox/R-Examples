hogg.test <-
function(x,y,...) {

	m<-length(x)
	n<-length(y)
	N<-m+n

	z<-c(x,y)
	R<-rank(z)

	q1<-Q1(z)
	q2<-Q2(z)

	U<-R/(N+1)

	if( q2 > 7 ) {
		a<-sign(U-0.5)
		scores<-'sign'
	} else {
		if( q1 > 2 ) {
			a<-ifelse(U<=0.5, 4*U-1.5,0.5)
			scores<-'bent'	
		} else if( q2 < 2 ) {
			a<-ifelse( U <= 0.25, 4*U-1, ifelse( U >= 0.75, 4*U-3, 0) )
			scores<-'light'
		} else {
			a<-sqrt(12)*(R/(N+1)-0.5)
			scores<-'Wilcoxon'
		}
	}

	W<-sum(a[1:m])

	mu<-m*mean(a)
	sigma<-sqrt(n*m/N*var(a))
	Z<-(W-mu)/sigma
	pval<-2*(1-pnorm(abs(Z)))

res<-list(statistic=W, p.value=pval,scores=scores)
class(res)<-'hogg.test'
res

}
