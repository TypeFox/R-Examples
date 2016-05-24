prop_mod <-
function(curr.index,data,maximal.mod,null.move.prob=0.5){

stay<-runif(1)

if(stay>null.move.prob){

big.X<-maximal.mod$x

full.labels<-attr(summary(maximal.mod)$terms,"term.labels")
full.order<-attr(summary(maximal.mod)$terms,"order")
full.factors<-attr(summary(maximal.mod)$terms,"factors")
full.terms<-attributes(big.X)$assign

dr<-drop_term(curr.index,data,maximal.mod)
ad<-add_term(curr.index,data,maximal.mod)

terms<-c(dr,ad)
types<-rep(c(1,2),c(length(dr),length(ad)))
total.choices<-length(terms)

if(length(terms)>0){

choose<-sample(x=1:total.choices,size=1)

termo<-terms[choose]
typo<-types[choose]

new.index<-curr.index
if(typo==1){ #drop move
new.index[full.terms==(1:length(full.labels))[full.labels==termo]]<-0
result<-list(new.index=new.index,type="drop",total.choices=total.choices,null.move.prob=null.move.prob)}

if(typo==2){ #add move
new.index[full.terms==(1:length(full.labels))[full.labels==termo]]<-1
result<-list(new.index=new.index,type="add",total.choices=total.choices,null.move.prob=null.move.prob)}} else{

result<-list(new.index=curr.index,type="null",total.choices=0,null.move.prob=null.move.prob)}} else{

result<-list(new.index=curr.index,type="null",total.choices=0,null.move.prob=null.move.prob)}

result}
