test.model.MSAR = function(data,simu,lag=NULL,id=1,u=NULL){
	
	N.samples = dim(data)[2]
	N.s = N.samples
	N.sim = dim(simu)[2]
	Bsim = N.sim/N.samples
	# Stationary distribution
	StaDist = NULL
	qqp = qqplot(data[,,id],simu[,,id],plot.it = FALSE)
    qq = trapz(qqp$x,abs(qqp$x-qqp$y))
    qq.sim = NULL
	for (k in 1:Bsim) {
		qqp.sim = qqplot(simu[,((k-1)*N.s+1):(k*N.s),id],simu[,,id],plot.it = FALSE)
		qq.sim[k] = trapz(qqp.sim$x,abs(qqp.sim$x-qqp$y))
	}
	StaDist$dd = qq
	StaDist$q.dd = quantile(qq.sim,probs=c(.05,.95))
	StaDist$p.value = 1-rank(c(qq,qq.sim))[1]/(Bsim+1)
	
	if (is.null(u)) {u=seq(min(data),max(data),length.out=30)}
	F = NULL
	F.sim = NULL
	F.s = matrix(0,length(u),Bsim)
	for (iu in 1:length(u)) {
		F[iu] = sum(data[,,id]<=u[iu])
		F.sim[iu] = sum(simu[,,id]<=u[iu])
		for (k in 1:Bsim) {F.s[iu,k] = sum(simu[,((k-1)*N.s+1):(k*N.s),id]<=u[iu])}
	}
	F = (F+1/2)/(length(data)+1)
	F.sim = (F.sim+1/2)/(length(simu)+1)
	F.s = (F.s+1/2)/(length(data)+1)
	S = trapz(u,(F-F.sim)^2/(F.sim*(1-F.sim)))
	S.s = NULL
	for ( k in 1:Bsim){S.s[k] = trapz(u,(F.sim-F.s[,k])^2/(F.sim*(1-F.sim)))}
	AD = list()
	AD$S = S
	AD$q.S = quantile(S.s,probs=c(.05,.95))
	AD$p.value =  1-rank(c(S,S.s))[1]/(Bsim+1)

	# Correlation function
	N.s = dim(data)[2]
	N.sim = dim(simu)[2]
	N = floor(N.sim/N.s)
	if (is.null(lag)) {lag = floor(dim(data)[1]/2)}
	C.data=matrix(0,lag,1)
	for (ex in 1:N.samples) {
		C = acf(data[,ex,id], lag.max=lag-1,type = "correlation",plot=FALSE)
		C.data=C.data+c(C$acf)
	}
	C.data = C.data/N.s
	C.sim=matrix(0,lag,N)
	for (k in 1:N) {
		for (ex in ((k-1)*N.s+1):(k*N.s)) {
			C = acf(simu[,ex,id], lag.max=lag-1, type = "correlation",plot=FALSE)
			C.sim[,k]=C.sim[,k]+c(C$acf)
		}
		C.sim[,k] = C.sim[,k]/N.s
	}
	Cor = list()
	dd = trapz(0:(lag-1),abs(C.data-apply(C.sim,1,mean)))
	dd.sim = NULL
	for (b in 1:Bsim) {
		dd.sim[b] = trapz(0:(lag-1),abs(C.sim[,b]-apply(C.sim,1,mean)))
	}
	Cor$dd = dd
	Cor$q.dd = quantile(dd.sim,probs=c(.05,.95))
	Cor$p.value = 1-rank(c(dd,dd.sim))[1]/(Bsim+1)	
	
 	# Cor1 = list()
 	# dd = sum(C.data^2)-sum(apply(C.sim,1,mean)^2)
 	# dd.sim = NULL
 	# for (b in 1:Bsim) {
 		# dd.sim[b] = sum(C.sim[,b]^2)-sum(apply(C.sim,1,mean)^2)
 	# }
 	# Cor1$dd = dd
 	# Cor1$q.dd = quantile(dd.sim,probs=c(.05,.95))
 	# Cor1$p.value = 1-rank(c(dd,dd.sim))[1]/(Bsim+1)

	
	# Up crossings
	u = seq(min(data),max(data),by=.3)
	Nu = matrix(0,length(u),N.samples)
	Nu.sim = matrix(0,length(u),N.sim)
	for (iu in 1:length(u)) {
		for (ex in 1:N.samples) {
			Nu[iu,ex] = sum(data[1:(T-1),ex,1]>u[iu] & data[2:T,ex,1]<u[iu] )
		}
		for (ex in 1:N.sim) {
			Nu.sim[iu,ex] = sum(simu[1:(T-1),ex,1]>u[iu] & simu[2:T,ex,1]<u[iu] )
		}	
	}
	Nu1 = apply(Nu,1,sum)/N.samples
	Nu1.sim = apply(Nu.sim,1,sum)/N.sim
	dNu = trapz(u,abs(Nu1-Nu1.sim))
	dNu.sim = NULL
	for (k in 1:Bsim) {
		Nu = apply(Nu.sim[,((k-1)*N.s+1):(k*N.s)],1,sum)/N.samples
		dNu.sim[k] = trapz(u,abs(Nu-Nu1.sim))
	}
    ENu = list()
	ENu$dd = dNu
	ENu$q.dd = quantile(dNu.sim,probs=c(.05,.95))
	r = rank(c(dNu,dNu.sim))[1]
	ENu$p.value = 1-r/(Bsim+1)
	
    return(list(StaDist =StaDist ,Cor=Cor,ENu=ENu,AD=AD))
}



