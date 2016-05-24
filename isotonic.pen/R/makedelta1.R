makedelta1 <-
function(x){
	n=length(x)
	rng=max(x)-min(x)
	mx=min(x)
	x=(x-mx)/rng
	x=round(x,8)
	x=x*rng+mx
	xu=sort(unique(x))
	nu=length(xu)
	delta=matrix(nrow=nu-1,ncol=n)
	for(i in 1:(nu-1)){
		delta[i,x<=xu[i]]=0;delta[i,x>xu[i]]=1
		delta[i,]=delta[i,]-mean(delta[i,])
	}
	delta
}
