rlba <-
function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
	n.with.extras=ceiling(n*(1+3*prod(pnorm(-vs))))
	drifts=matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE)
	if (truncdrifts) {
		repeat {
			drifts=rbind(drifts,matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE))
			tmp=apply(drifts,1,function(x) any(x>0))
			drifts=drifts[tmp,]
			if (nrow(drifts)>=n) break
		}
	}
	drifts=drifts[1:n,]
	drifts[drifts<0]=0
	starts=matrix(runif(min=0,max=A,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
	ttf=t((b-t(starts)))/drifts
	rt=apply(ttf,1,min)+t0+runif(min=-st0/2,max=+st0/2,n=n)
	resp=apply(ttf,1,which.min)
	data.frame(resp=resp,rt=rt)
}
