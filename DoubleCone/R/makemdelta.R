makemdelta <-
function(x){
	n=length(x)
	xs=sort(x)
	xu=1:n*0
# find unique x values
	nu=1;xu[1]=xs[1];sm=1e-12
	for(i in 2:n){
		if(xs[i]>xs[i-1]+sm){
			nu=nu+1
			xu[nu]=xs[i]
		}
	}
	delta=matrix(nrow=nu-1,ncol=n)
	for(i in 1:(nu-1)){
		delta[i,x<=xu[i]]=0;delta[i,x>xu[i]]=1
		delta[i,]=delta[i,]-mean(delta[i,])
	}
	delta
}
