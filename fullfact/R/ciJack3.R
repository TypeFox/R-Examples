ciJack3 <-
function(comp,full,remain=NULL,level=95,rnd_r=3,rnd_p=1,trait=NULL) {
  if (missing(full)) stop("Need the observed values using full data")
  cia<- (100-level)/100/2
  n<- nrow(comp)
  mater<- grep("maternal", colnames(comp))
  add<- grep("additive", colnames(comp))
  nonadd<- grep("nonadd", colnames(comp))
if (!is.null(remain)) {
  play<- matrix(0,ncol=1,nrow=length(remain))  #remaining columns
   for (i in 1:length(remain)) { play[i,]<- grep(paste(remain[i]), colnames(comp)) } }
  perc<- matrix(0,ncol=length(remain)+3,nrow=nrow(comp))
  perc[,1]<- 100*comp[,add]/comp$Total
  perc[,2]<- 100*comp[,nonadd]/comp$Total
  perc[,3]<- 100*comp[,mater]/comp$Total
if (!is.null(remain)) {
  for (i in 1:length(remain)) { perc[,(i+3)]<- 100*comp[,play[i,]]/comp$Total } }
if (is.null(remain)) {
 full_r<- do.call("rbind", replicate(n,full[1:3],simplify=F)) #raw
 full_p1<- matrix(0,ncol=3,nrow=1) }  #percentage
if (!is.null(remain)) {
 full_r<- do.call("rbind", replicate(n,full[-4],simplify=F)) #no total
 full_p1<- matrix(0,ncol=3+length(remain),nrow=1) }
  full_p1[,1]<- 100*full[1]/full[4] #additive
  full_p1[,2]<- 100*full[2]/full[4] #non-additive
  full_p1[,3]<- 100*full[3]/full[4] #maternal
if (!is.null(remain)) {
  for (i in 1:length(remain)) { full_p1[,(3+i)]<- 100*full[(4+i)]/full[4] } }
full_p<- do.call("rbind", replicate(n,full_p1,simplify=F))
if (is.null(remain)) {
  pseudo_r<- n*full_r - (n-1)*comp[,c(add,nonadd,mater)]
  pseudo_p<- n*full_p - (n-1)*perc[,c(1:3)] }
if (!is.null(remain)) {
  pseudo_r<- n*full_r - (n-1)*comp[,c(add,nonadd,mater,play)]
  pseudo_p<- n*full_p - (n-1)*perc[,c(1:(3+length(remain)))] }
lwr_r<- apply(pseudo_r, 2, mean)- qt(1-cia,n-1)*sqrt(apply(pseudo_r, 2, var)/n)
med_r<- apply(pseudo_r, 2, mean)
upp_r<- apply(pseudo_r, 2, mean) + qt(1-cia,n-1)*sqrt(apply(pseudo_r, 2, var)/n)
lwr_p<- apply(pseudo_p, 2, mean)- qt(1-cia,n-1)*sqrt(apply(pseudo_p, 2, var)/n)
med_p<- apply(pseudo_p, 2, mean)
upp_p<- apply(pseudo_p, 2, mean) + qt(1-cia,n-1)*sqrt(apply(pseudo_p, 2, var)/n)
if (is.null(remain)) { ci<- matrix(0,ncol=4,nrow=3); ci_p<- matrix(0,ncol=4,nrow=3) }
if (!is.null(remain)) { ci<- matrix(0,ncol=4,nrow=3+length(remain)); ci_p<- matrix(0,ncol=4,nrow=3+length(remain)) }
col_names1<- c("component","lower","median","upper") #know column names
ci[,1][1:3]<- c("additive","nonadd","maternal") #known labels
ci_p[,1][1:3]<- c("additive","nonadd","maternal")
ci[,2][1:3]<- lwr_r[1:3]; ci_p[,2][1:3]<- lwr_p[1:3]  #known lower
ci[,3][1:3]<- med_r[1:3]; ci_p[,3][1:3]<- med_p[1:3] #known mean
ci[,4][1:3]<- upp_r[1:3]; ci_p[,4][1:3]<- upp_p[1:3]   #known upper
if (!is.null(remain)) { for (i in 1:length(remain)) {
  ci[,1][(3+i)]<- paste(remain[i]); ci_p[,1][(3+i)]<- paste(remain[i])
  ci[(3+i),][2:4]<- c(lwr_r[3+i],med_r[3+i],upp_r[3+i]); ci_p[(3+i),][2:4]<- c(lwr_p[3+i],med_p[3+i],upp_p[3+i]) }  }
  ci[,2:4]<- round(as.numeric(ci[,2:4]),rnd_r); ci_p[,2:4]<- round(as.numeric(ci_p[,2:4]),rnd_p)
  ci<- as.data.frame(ci); ci_p<- as.data.frame(ci_p)
  colnames(ci)<- c("component","lower","mean","upper"); colnames(ci_p)<- c("component","lower","mean","upper")
  if (is.null(trait)) { ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
  if (!is.null(trait)) { ci$trait<- as.factor(trait); ci_p$trait<- as.factor(trait);
    ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
}
