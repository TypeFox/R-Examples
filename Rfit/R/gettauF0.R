gettauF0<-function(ehat,p=0,scores=Rfit::wscores,delta=0.8,hparm=2,...) {

n<-length(ehat)

sciint<-mad(ehat)

	.Fortran('nscale',
		as.integer(n),
		as.double(.Machine$double.eps^0.25),
		as.double(delta),
		as.double(sciint),
		as.integer(0),
		as.integer(p),
		as.double(ehat),
		as.integer(order(ehat)),
		as.double(getScores(scores,seq_len(n)/(n+1))),
		as.double(rep.int(0,n)),
		as.double(rep.int(0,n)),
		tauhat=as.double(0),
		as.double(rep.int(0,5)),
		as.integer(0),
		as.integer(1000),
		as.double(0),
		as.double(getScoresDeriv(scores,seq_len(n)/(n+1))),
		as.double(hparm),
		PACKAGE='Rfit'
	)$tauhat
		
		
}
