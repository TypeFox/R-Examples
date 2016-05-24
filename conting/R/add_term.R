add_term <-
function(curr.index,data,maximal.mod){

big.X<-maximal.mod$x
full.terms<-attributes(big.X)$assign

full.labels<-attr(summary(maximal.mod)$terms,"term.labels")
full.order<-attr(summary(maximal.mod)$terms,"order")
full.factors<-attr(summary(maximal.mod)$terms,"factors")

uni<-unique(full.terms[curr.index==1])
uni<-uni[uni>0]
curr.labels<-full.labels[uni]
curr.order<-full.order[uni]
curr.factors<-full.factors[,uni]

K<-length(full.labels[full.order==1])

can_add<-c()
if(!all(curr.index==1)){
pot.labels<-c()
pot.order<-c()
pot<-c()
for(ttt in 1:length(full.labels)){
if(!any(full.labels[ttt]==curr.labels)){
pot.labels<-c(pot.labels,full.labels[ttt])
pot.order<-c(pot.order,full.order[ttt])
pot<-c(pot,ttt)}}

pot.factors<-matrix(full.factors[-1,pot],nrow=K)

can_add<-pot.labels[pot.order==2]
candos<-(1:length(pot.labels))[pot.order>2]
candos<-candos[pot.order[candos]<=(max(curr.order)+1)]
for(ttt in candos){
combos<-combinations(n=pot.order[ttt],r=pot.order[ttt]-1,v=(1:K)[pot.factors[,ttt]==1])
ok<-c()
for(q in 1:dim(combos)[1]){
int<-rep(0,K)
int[combos[q,]]<-1
run<-as.numeric(apply(matrix(rep(int,dim(curr.factors)[2]),nrow=K,byrow=FALSE)==curr.factors[-1,],2,all))
ok[q]<-ifelse(sum(run)==0,0,1)}
if(all(ok==1)){
can_add<-c(can_add,pot.labels[ttt])}}}

can_add}
