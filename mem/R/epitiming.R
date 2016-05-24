epitiming <-
function(i.data,i.n.values=5,i.method=2,i.param=2.8){
	loc.epi<-localizar.epidemia(i.data,i.n.values,i.method,i.param)
  names(loc.epi)<-c("i.data","map.curve","optimum.map","pre.epi","post.epi")
  loc.epi$call<-match.call()
	class(loc.epi)<-"epidemic"
	return(loc.epi)
}
