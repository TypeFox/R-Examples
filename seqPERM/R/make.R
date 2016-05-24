make <-
function(r1,r2,k,v){
S<-seq(r1,r2,1)
L<-length(S)
	OOO<-function(L,k,v){
		OO<-function(L,k,v){
			s1<-lapply(1:L^(v-1),function(i) c(S[k]))
			ss1<-do.call("c",s1)
		return(ss1)}
	SS<-lapply(1:L, function(i) c(OO(L,i,k)))
	SSS<-do.call("c",SS)
	return(SSS)}
l<-lapply(1:L^(v-k),function(i) c(OOO(L,k,v)))
ll<-do.call("c",l)
return(ll)}
