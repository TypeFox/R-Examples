drop_term <-
function(curr.index,data,maximal.mod){

big.X<-maximal.mod$x							## maximal design matrix
full.terms<-attributes(big.X)$assign					## terms in big.X

uni<-unique(full.terms[curr.index==1])					## terms in current model
uni<-uni[uni>0]								## non-interecpt terms in current model
term.labels<-attr(summary(maximal.mod)$terms,"term.labels")[uni]		## term labels in current model
term.order<-attr(summary(maximal.mod)$terms,"order")[uni]			## order of terms in current model
term.factors<-attr(summary(maximal.mod)$terms,"factors")[,uni]		## constituent main effect terms of terms in 									## current model

K<-length(term.labels[term.order==1])					## number of main effect terms
can_drop<-c()								## vector containing terms we can drop
if(max(term.order)>1){							## if we are not in independence model
can_drop<-term.labels[term.order==max(term.order)]			## we can drop all highest order terms
candos<-(1:length(term.labels[term.order<max(term.order)]))[-(1:K)]	## all non highest order terms and non main effect terms are 		
for(ttt in candos){							## candidates to drop	
tls<-(1:K)[term.factors[-1,ttt]==1]					## terms that depend on this term
ntls<-(1:K)[term.factors[-1,ttt]==0]					## terms that do not depend on this term
ok<-c()
for(q in 1:length(ntls)){
TLS<-c(tls,ntls[q])
int<-rep(0,K)
int[TLS]<-1
run<-as.numeric(apply(matrix(rep(int,dim(term.factors)[2]),nrow=K,byrow=FALSE)==term.factors[-1,],2,all)) 
ok[q]<-ifelse(sum(run)==0,1,0)}
if(all(ok==1)){
can_drop<-c(can_drop,term.labels[ttt])}}}

can_drop}
