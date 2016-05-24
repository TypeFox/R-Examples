med1way.crit <-
function(n,alpha=.05,iter=1000,TEST=NA,SEED=TRUE){
#
#  Determine the critical value for the function
#  med1way, assuming normality, based on the sample
#  sizes in n.
#
J<-length(n)
x<-list()
w<-vector("numeric",J)
xbar<-vector("numeric",J)
if(SEED)set.seed(2)
chk<-NA
grp<-c(1:J)
for (it in 1:iter){
for(j in 1:J){
x[[j]]<-rnorm(n[j])
w[j]<-1/msmedse(x[[grp[j]]])^2
xbar[j]<-median(x[[grp[j]]])
n[j]<-length(x[[grp[j]]])
}
u<-sum(w)
xtil<-sum(w*xbar)/u
chk[it]<-sum(w*(xbar-xtil)^2)/(J-1)
}
chk<-sort(chk)
iv<-round((1-alpha)*iter)
crit.val<-chk[iv]
pval<-NA
if(!is.na(TEST))pval<-sum((TEST<=chk))/iter
list(crit.val=crit.val,p.value=pval)
}
