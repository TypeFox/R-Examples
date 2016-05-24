rank.test<-function(x,y,alternative="two.sided",scores=Rfit::wscores,conf.int=FALSE,conf.level=0.95) {
     z = c(x,y)
     n1 = length(x)
     n2 = length(y)
     n = n1 + n2
     asc = getScores(scores, 1:n/(n+1))
     rs = rank(z)/(n+1)
     asg = getScores(scores,rs)
     Sphi = sum(asg[(n1+1):n])
     sigphi = sqrt(((n1*n2)/(n*(n-1)))*sum(asc^2))
     zphi = Sphi/sigphi
	pvalue<-switch(alternative,
			two.sided = 2*(1 - pnorm(abs(zphi))),
			less = pnorm(zphi),
			greater = 1 - pnorm(zphi)
		)
     res<-list(Sphi=Sphi,statistic=zphi,p.value=pvalue)
	if( conf.int ) {
		w<-c(rep(0,n1),rep(1,n2))
		fit<-rfit(z~w)
		estse<-coef(summary(fit))[2,1:2]
		alpha<-1-conf.level

		tcvs<-switch(alternative,
			two.sided = qt(1-alpha/2,n-2)*c(-1,1),
			less = c(-Inf,qt(1-alpha,n-2)),
			greater = c(qt(alpha,n-2),Inf)
		)

		res$conf.int<-estse[1]+tcvs*estse[2]
		res$estimate<-estse[1]

	}
	class(res)<-'rank.test'
	res
}
