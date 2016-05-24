ciJack2 <-
function(comp,full,level=95,rnd_r=3,rnd_p=1,position=NULL,block=NULL,trait=NULL) {
  if (missing(full)) stop("Need the observed values using full data")
  cia<- (100-level)/100/2
  n<- nrow(comp)
  mater<- grep("maternal", colnames(comp))
  add<- grep("additive", colnames(comp))
  nonadd<- grep("nonadd", colnames(comp))
  if (!is.null(position)) { pos<- grep(paste(position), colnames(comp)) }
  if (!is.null(block)) { bloc<- grep(paste(block), colnames(comp)) }
  comp$p_mat<- 100*comp[,mater]/comp$Total
  comp$p_add<- 100*comp[,add]/comp$Total
  comp$p_na<- 100*comp[,nonadd]/comp$Total
  if (!is.null(position)) { comp$p_pos<- 100*comp[,pos]/comp$Total }
  if (!is.null(block)) { comp$p_bloc<- 100*comp[,bloc]/comp$Total }
if (is.null(position) && is.null(block)) {
 full_r<- do.call("rbind", replicate(n,full[1:3],simplify=F)) #raw
 full_p1<- data.frame(additive=100*full[1]/full[4], nonadd=100*full[2]/full[4], maternal=100*full[3]/full[4]) } #percentage
if (!is.null(position) && is.null(block)) {
 full_r<- do.call("rbind", replicate(n,full[c(1:3,5)],simplify=F))
 full_p1<- data.frame(additive=100*full[1]/full[4],nonadd=100*full[2]/full[4],maternal=100*full[3]/full[4],
   position=100*full[5]/full[4]) }
if (is.null(position) && !is.null(block)) {
 full_r<- do.call("rbind", replicate(n,full[c(1:3,5)],simplify=F))
 full_p1<- data.frame(additive=100*full[1]/full[4],nonadd=100*full[2]/full[4],maternal=100*full[3]/full[4],
   block=100*full[5]/full[4]) }
if (!is.null(position) && !is.null(block)) {
 full_r<- do.call("rbind", replicate(n,full[c(1:3,5,6)],simplify=F))
 full_p1<- data.frame(additive=100*full[1]/full[4],nonadd=100*full[2]/full[4],maternal=100*full[3]/full[4],
   position=100*full[5]/full[4],block=100*full[6]/full[4])  }
full_p<- do.call("rbind", replicate(n,full_p1,simplify=F))
if (is.null(position) && is.null(block)) {
  pseudo_r<- n*full_r - (n-1)*comp[,c(add,nonadd,mater)]
  pseudo_p<- n*full_p - (n-1)*cbind(comp$p_add,comp$p_na,comp$p_mat) }
if (!is.null(position) && is.null(block)) {
  pseudo_r<- n*full_r - (n-1)*comp[,c(add,nonadd,mater,pos)]
  pseudo_p<- n*full_p - (n-1)*cbind(comp$p_add,comp$p_na,comp$p_mat,comp$p_pos) }
if (is.null(position) && !is.null(block)) {
  pseudo_r<- n*full_r - (n-1)*comp[,c(add,nonadd,mater,bloc)]
  pseudo_p<- n*full_p - (n-1)*cbind(comp$p_add,comp$p_na,comp$p_mat,comp$p_bloc) }
if (!is.null(position) && !is.null(block)) {
  pseudo_r<- n*full_r - (n-1)*comp[,c(add,nonadd,mater,pos,bloc)]
  pseudo_p<- n*full_p - (n-1)*cbind(comp$p_add,comp$p_na,comp$p_mat,comp$p_pos,comp$p_bloc) }
lwr_r<- colMeans(pseudo_r)- qt(1-cia,n-1)*sqrt(apply(pseudo_r, 2, var)/n)
upp_r<- colMeans(pseudo_r) + qt(1-cia,n-1)*sqrt(apply(pseudo_r, 2, var)/n)
lwr_p<- colMeans(pseudo_p)- qt(1-cia,n-1)*sqrt(apply(pseudo_p, 2, var)/n)
upp_p<- colMeans(pseudo_p) + qt(1-cia,n-1)*sqrt(apply(pseudo_p, 2, var)/n)
if (is.null(position) && is.null(block)) {
  ci<- data.frame(component=c("additive","nonadd","maternal"),lower= lwr_r, mean=colMeans(pseudo_r), upper= upp_r)
  rownames(ci)<- c(1,2,3)
  ci_p<- data.frame(component=c("additive","nonadd","maternal"),lower= lwr_p, mean=colMeans(pseudo_p), upper= upp_p)
  rownames(ci_p)<- c(1,2,3) }
if (!is.null(position) && is.null(block)) {
  ci<- data.frame(component=c("additive","nonadd","maternal",paste(position)),lower= lwr_r, mean=colMeans(pseudo_r), upper= upp_r)
  rownames(ci)<- c(1,2,3,4)
  ci_p<- data.frame(component=c("additive","nonadd","maternal",paste(position)),lower= lwr_p, mean=colMeans(pseudo_p), upper= upp_p)
  rownames(ci_p)<- c(1,2,3,4) }
if (is.null(position) && !is.null(block)) {
  ci<- data.frame(component=c("additive","nonadd","maternal",paste(block)),lower= lwr_r, mean=colMeans(pseudo_r), upper= upp_r)
  rownames(ci)<- c(1,2,3,4)
  ci_p<- data.frame(component=c("additive","nonadd","maternal",paste(block)),lower= lwr_p, mean=colMeans(pseudo_p), upper= upp_p)
  rownames(ci_p)<- c(1,2,3,4) }
if (!is.null(position) && !is.null(block)) {
  ci<- data.frame(component=c("additive","nonadd","maternal",paste(position),paste(block)),lower= lwr_r,
   mean=colMeans(pseudo_r), upper= upp_r); rownames(ci)<- c(1,2,3,4,5)
  ci_p<- data.frame(component=c("additive","nonadd","maternal",paste(position),paste(block)),lower= lwr_p,
   mean=colMeans(pseudo_p), upper= upp_p);rownames(ci_p)<- c(1,2,3,4,5) }
ci[,2:4]<- round(ci[,2:4],rnd_r)
ci_p[,2:4]<- round(ci_p[,2:4],rnd_p)
  if (is.null(trait)) { ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
  if (!is.null(trait)) { ci$trait<- as.factor(trait); ci_p$trait<- as.factor(trait);
    ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
}
