ciMANA3 <-
function(comp,level=95,rnd_r=3,rnd_p=1,bias=NULL,accel=NULL,remain=NULL,trait=NULL) {
  cia<- (100-level)/100/2
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
if (!is.null(bias)) {
  z0_mat <- qnorm(mean(comp[,mater] < bias[3]))
  z0_add <- qnorm(mean(comp[,add] < bias[1]))
  z0_na <- qnorm(mean(comp[,nonadd] < bias[2]))
if (!is.null(remain)) {
  z0_play<- matrix(0,ncol=1,nrow=length(remain))  #remaining z0
   for (i in 1:length(remain)) { z0_play[i,]<- qnorm(mean(comp[,play[i,]] < bias[3+i])) } }
} #end bias play
if (is.null(accel)) { a_mat<- 0; a_add<- 0; a_na<- 0
  if (!is.null(remain)) { a_play<- matrix(0,ncol=1,nrow=length(remain)) } }  #remaining a, all zero
if (!is.null(accel)) {
  mater2<- grep("maternal", colnames(accel))
  add2<- grep("additive", colnames(accel))
  nonadd2<- grep("nonadd", colnames(accel))
if (!is.null(remain)) {
  play2<- matrix(0,ncol=1,nrow=length(remain))  #remaining columns
   for (i in 1:length(remain)) { play2[i,]<- grep(paste(remain[i]), colnames(accel)) } }
if (!is.null(remain)) {
  for (i in 1:length(remain)) { accel[,play2[i,]]<- 100*accel[,play2[i,]]/accel$Total } }
  a_mat <- sum((bias[3]-accel[,mater2])^3)/(6*sum((bias[3]-accel[,mater2])^2)^(3/2))
  a_add <- sum((bias[1]-accel[,add2])^3)/(6*sum((bias[1]-accel[,add2])^2)^(3/2))
  a_na <- sum((bias[2]-accel[,nonadd2])^3)/(6*sum((bias[2]-accel[,nonadd2])^2)^(3/2))
if (!is.null(remain)) { a_play<- matrix(0,ncol=1,nrow=length(remain))
for (i in 1:length(remain)) {
  a_play[i,]<- sum((bias[3+i]-accel[,play[i,]])^3)/(6*sum((bias[3+i]-accel[,play[i,]])^2)^(3/2)) }
  } #end remain
} #end acceleration
if (is.null(remain)) { ci<- matrix(0,ncol=4,nrow=3); ci_p<- matrix(0,ncol=4,nrow=3) }
if (!is.null(remain)) { ci<- matrix(0,ncol=4,nrow=3+length(remain)); ci_p<- matrix(0,ncol=4,nrow=3+length(remain)) }
  col_names1<- c("component","lower","median","upper") #know column names
  ci[,1][1:3]<- c("additive","nonadd","maternal") #known labels
  ci_p[,1][1:3]<- c("additive","nonadd","maternal")
  ci[,2][1:3]<- c(quantile(comp[,add],cia),quantile(comp[,nonadd],cia),quantile(comp[,mater],cia)) #known lower
  ci[,3][1:3]<- c(quantile(comp[,add],0.5),quantile(comp[,nonadd],0.5),quantile(comp[,mater],0.5))  #known median
  ci[,4][1:3]<- c(quantile(comp[,add],1-cia),quantile(comp[,nonadd],1-cia),quantile(comp[,mater],1-cia)) #known upper
if (!is.null(remain)) { for (i in 1:length(remain)) {
  ci[,1][(3+i)]<- paste(remain[i])
  ci[(3+i),][2:4]<- c(quantile(comp[,play[i,]],cia),quantile(comp[,play[i,]],0.5),quantile(comp[,play[i,]],1-cia)) }  }
  ci_p[,2][1:3]<- c(quantile(perc[,1],cia),quantile(perc[,2],cia),quantile(perc[,3],cia)) #known lower
  ci_p[,3][1:3]<- c(quantile(perc[,1],0.5),quantile(perc[,2],0.5),quantile(perc[,3],0.5))  #known median
  ci_p[,4][1:3]<- c(quantile(perc[,1],1-cia),quantile(perc[,2],1-cia),quantile(perc[,3],1-cia)) #known upper
if (!is.null(remain)) { for (i in 1:length(remain)) {
  ci_p[,1][(3+i)]<- paste(remain[i])
  ci_p[(3+i),][2:4]<- c(quantile(perc[,(3+i)],cia),quantile(perc[,(3+i)],0.5),quantile(perc[,(3+i)],1-cia)) }  }
if (!is.null(bias)) {
ql_mat <- pnorm(z0_mat+(z0_mat+qnorm(cia))/(1-a_mat*(z0_mat+qnorm(cia))))
ql_add <- pnorm(z0_add+(z0_add+qnorm(cia))/(1-a_add*(z0_add+qnorm(cia))))
ql_na <- pnorm(z0_na+(z0_na+qnorm(cia))/(1-a_na*(z0_na+qnorm(cia))))
md_mat <- pnorm(z0_mat+(z0_mat+qnorm(0.50))/(1-a_mat*(z0_mat+qnorm(0.50))))
md_add <- pnorm(z0_add+(z0_add+qnorm(0.50))/(1-a_add*(z0_add+qnorm(0.50))))
md_na <- pnorm(z0_na+(z0_na+qnorm(0.50))/(1-a_na*(z0_na+qnorm(0.50))))
qu_mat <- pnorm(z0_mat+(z0_mat+qnorm(1-cia))/(1-a_mat*(z0_mat+qnorm(1-cia))))
qu_add <- pnorm(z0_add+(z0_add+qnorm(1-cia))/(1-a_add*(z0_add+qnorm(1-cia))))
qu_na <- pnorm(z0_na+(z0_na+qnorm(1-cia))/(1-a_na*(z0_na+qnorm(1-cia))))
if (!is.null(remain)) {
  ql_play<- matrix(0,ncol=1,nrow=length(remain))
  md_play<- matrix(0,ncol=1,nrow=length(remain))
  qu_play<- matrix(0,ncol=1,nrow=length(remain))
for (i in 1:length(remain)) {
  ql_play[i,]<- pnorm(z0_play[i,]+(z0_play[i,]+qnorm(cia))/(1-a_play[i,]*(z0_play[i,]+qnorm(cia))))
  md_play[i,]<- pnorm(z0_play[i,]+(z0_play[i,]+qnorm(0.5))/(1-a_play[i,]*(z0_play[i,]+qnorm(0.5))))
  qu_play[i,]<- pnorm(z0_play[i,]+(z0_play[i,]+qnorm(1-cia))/(1-a_play[i,]*(z0_play[i,]+qnorm(1-cia))))
  }  #end remain loop
} #end bias constants
if (is.null(remain)) { ci2<- matrix(0,ncol=4,nrow=3); ci2_p<- matrix(0,ncol=4,nrow=3)  }
if (!is.null(remain)) { ci2<- matrix(0,ncol=4,nrow=3+length(remain)); ci2_p<- matrix(0,ncol=4,nrow=3+length(remain)) }
  col_names<- c("component","lower","median","upper") #known column names
  ci2[,1][1:3]<- c("additive","nonadd","maternal") #known labels
  ci2_p[,1][1:3]<- c("additive","nonadd","maternal")
  ci2[,2][1:3]<- c(quantile(comp[,add],ql_add),quantile(comp[,nonadd],ql_na),quantile(comp[,mater],ql_mat)) #known lower
  ci2[,3][1:3]<- c(quantile(comp[,add],md_add),quantile(comp[,nonadd],md_na),quantile(comp[,mater],md_mat)) #known median
  ci2[,4][1:3]<- c(quantile(comp[,add],qu_add),quantile(comp[,nonadd],qu_na),quantile(comp[,mater],qu_mat)) #known upper
if (!is.null(remain)) { for (i in 1:length(remain)) {
  ci2[,1][(3+i)]<- paste(remain[i])
  ci2[(3+i),][2:4]<- c(quantile(comp[,play[i,]],ql_play[i,]),quantile(comp[,play[i,]],md_play[i,]),quantile(comp[,play[i,]],qu_play[i,]))
  } #end loop
} #end remain
  ci2_p[,2][1:3]<- c(quantile(perc[,1],ql_add),quantile(perc[,2],ql_na),quantile(perc[,3],ql_mat)) #known lower
  ci2_p[,3][1:3]<- c(quantile(perc[,1],md_add),quantile(perc[,2],md_na),quantile(perc[,3],md_mat)) #known median
  ci2_p[,4][1:3]<- c(quantile(perc[,1],qu_add),quantile(perc[,2],qu_na),quantile(perc[,3],qu_mat)) #known upper
if (!is.null(remain)) { for (i in 1:length(remain)) {
  ci2_p[,1][(3+i)]<- paste(remain[i])
  ci2_p[(3+i),][2:4]<- c(quantile(perc[,(3+i)],ql_play[i,]),quantile(perc[,(3+i)],md_play[i,]),quantile(perc[,(3+i)],qu_play[i,]))
  } #end loop
} #end remain
if (z0_add == Inf | z0_add == -Inf | z0_na == Inf | z0_na == -Inf | z0_mat == Inf | z0_mat == -Inf)
  { if (is.null(remain)) { ci2.1<- matrix(NA,ncol=5,nrow=3);ci2.1p<- matrix(NA,ncol=5,nrow=3) }
    if (!is.null(remain)) { ci2.1<- matrix(NA,ncol=5,nrow=3+length(remain)); ci2.1p<- matrix(NA,ncol=5,nrow=3+length(remain)) }
    ci2.1[,1:4]<- ci2[,1:4]; ci2<- ci2.1
    ci2.1p[,1:4]<- ci2_p[,1:4]; ci2_p<- ci2.1p
    col_names<- c("component","lower","median","upper","change") }
if (!is.null(remain)) { chg_test<- matrix(0,ncol=1,nrow=length(remain))
  for (i in 1:length(remain)) {
  if (z0_play[i,] == Inf | z0_play[i,] == -Inf) { chg_test[i,]<-1 }  }
 if (sum(chg_test) > 0 ) {
  ci2.1<- matrix(0,ncol=5,nrow=3+length(remain));ci2.1p<- matrix(0,ncol=5,nrow=3+length(remain))
  ci2.1[,1:4]<- ci2[,1:4]; ci2<- ci2.1
  ci2.1p[,1:4]<- ci2_p[,1:4]; ci2_p<- ci2.1p
  col_names<- c("component","lower","median","upper","change") } }
if (z0_add == Inf | z0_add == -Inf) { ci2[1,][2:4]<- ci[1,][2:4]; ci2[1,][5] <- "bias fail"
  ci2_p[1,][2:4]<- ci_p[1,][2:4]; ci2_p[1,][5] <- "bias fail" }
if (z0_na == Inf | z0_na == -Inf) { ci2[2,][2:4]<- ci[2,][2:4]; ci2[2,][5] <- "bias fail"
  ci2_p[2,][2:4]<- ci_p[2,][2:4]; ci2_p[2,][5] <- "bias fail"  }
if (z0_mat == Inf | z0_mat == -Inf) { ci2[3,][2:4]<- ci[3,][2:4]; ci2[3,][5]<- "bias fail"
  ci2_p[3,][2:4]<- ci_p[3,][2:4]; ci2_p[3,][5]<- "bias fail"   }
if (!is.null(remain) && sum(chg_test) > 0) { for (i in 1:length(remain)) {
   ci2[(3+i),][2:4]<- ci[(3+i),][2:4]; ci2[(3+i),][5]<- "bias fail"
   ci2_p[(3+i),][2:4]<- ci_p[(3+i),][2:4]; ci2_p[(3+i),][5]<- "bias fail" } }
  ci2[,2:4]<- round(as.numeric(ci2[,2:4]),rnd_r)
  ci2_p[,2:4]<- round(as.numeric(ci2_p[,2:4]),rnd_p)
  ci2<- as.data.frame(ci2); colnames(ci2)<- col_names
  ci2_p<- as.data.frame(ci2_p); colnames(ci2_p)<- col_names
} #end null bias
  ci[,2:4]<- round(as.numeric(ci[,2:4]),rnd_r)
  ci_p[,2:4]<- round(as.numeric(ci_p[,2:4]),rnd_p)
  ci<- as.data.frame(ci);colnames(ci)<- col_names1
  ci_p<- as.data.frame(ci_p);colnames(ci_p)<- col_names1
  if (is.null(trait) == T && is.null(bias)) { ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
  if (is.null(trait) == T && !is.null(bias)) { ci_obj<- list(raw=ci2,percentage=ci2_p); return(ci_obj) }
  if (is.null(trait) == F && is.null(bias)) { ci$trait<- as.character(trait); ci_p$trait<- as.character(trait)
   ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
  if (is.null(trait) == F && !is.null(bias)) { ci2$trait<- as.character(trait); ci2_p$trait<- as.character(trait)
   ci_obj<- list(raw=ci2,percentage=ci2_p); return(ci_obj) }
}
