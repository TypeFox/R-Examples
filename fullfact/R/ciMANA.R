ciMANA <-
function(comp,level=95,rnd_r=3,rnd_p=1,bias=NULL,accel=NULL,trait=NULL) {
  cia<- (100-level)/100/2
  mater<- grep("maternal", colnames(comp))
  add<- grep("additive", colnames(comp))
  nonadd<- grep("nonadd", colnames(comp))
  comp$p_mat<- 100*comp[,mater]/comp$Total
  comp$p_add<- 100*comp[,add]/comp$Total
  comp$p_na<- 100*comp[,nonadd]/comp$Total
if (!is.null(bias)) {
  z0_mat <- qnorm(mean(comp[,mater] < bias[3]))
  z0_add <- qnorm(mean(comp[,add] < bias[1]))
  z0_na <- qnorm(mean(comp[,nonadd] < bias[2])) }
if (is.null(accel)) { a_mat<- 0; a_add<- 0; a_na<- 0 }
if (!is.null(accel)) {
  mater2<- grep("maternal", colnames(accel))
  add2<- grep("additive", colnames(accel))
  nonadd2<- grep("nonadd", colnames(accel))
  a_mat <- sum((bias[3]-accel[,mater2])^3)/(6*sum((bias[3]-accel[,mater2])^2)^(3/2))
  a_add <- sum((bias[1]-accel[,add2])^3)/(6*sum((bias[1]-accel[,add2])^2)^(3/2))
  a_na <- sum((bias[2]-accel[,nonadd2])^3)/(6*sum((bias[2]-accel[,nonadd2])^2)^(3/2))
} #end acceleration
  ci<- data.frame(component=c("additive","nonadd","maternal"),
    lower=c(quantile(comp[,add],cia),quantile(comp[,nonadd],cia),quantile(comp[,mater],cia)),
    median= c(quantile(comp[,add],0.5),quantile(comp[,nonadd],0.5),quantile(comp[,mater],0.5)),
    upper= c(quantile(comp[,add],1-cia),quantile(comp[,nonadd],1-cia),quantile(comp[,mater],1-cia)))
  ci<- cbind(ci$component,round(ci[,2:4],rnd_r));colnames(ci)[1]<- "component"
  rownames(ci)<- 1:3
  ci_p<- data.frame(component=c("additive","nonadd","maternal"),
    lower=c(quantile(comp$p_add,cia),quantile(comp$p_na,cia),quantile(comp$p_mat,cia)),
    median= c(quantile(comp$p_add,0.5),quantile(comp$p_na,0.5),quantile(comp$p_mat,0.5)),
    upper= c(quantile(comp$p_add,1-cia),quantile(comp$p_na,1-cia),quantile(comp$p_mat,1-cia)))
  ci_p<- cbind(ci_p$component,round(ci_p[,2:4],rnd_p));colnames(ci_p)[1]<- "component"
  rownames(ci_p)<- 1:3
if (!is.null(bias)) { #CI (adjusted)
ql_mat <- pnorm(z0_mat+(z0_mat+qnorm(cia))/(1-a_mat*(z0_mat+qnorm(cia))))
ql_add <- pnorm(z0_add+(z0_add+qnorm(cia))/(1-a_add*(z0_add+qnorm(cia))))
ql_na <- pnorm(z0_na+(z0_na+qnorm(cia))/(1-a_na*(z0_na+qnorm(cia))))
md_mat <- pnorm(z0_mat+(z0_mat+qnorm(0.50))/(1-a_mat*(z0_mat+qnorm(0.50))))
md_add <- pnorm(z0_add+(z0_add+qnorm(0.50))/(1-a_add*(z0_add+qnorm(0.50))))
md_na <- pnorm(z0_na+(z0_na+qnorm(0.50))/(1-a_na*(z0_na+qnorm(0.50))))
qu_mat <- pnorm(z0_mat+(z0_mat+qnorm(1-cia))/(1-a_mat*(z0_mat+qnorm(1-cia))))
qu_add <- pnorm(z0_add+(z0_add+qnorm(1-cia))/(1-a_add*(z0_add+qnorm(1-cia))))
qu_na <- pnorm(z0_na+(z0_na+qnorm(1-cia))/(1-a_na*(z0_na+qnorm(1-cia))))
  ci2<- data.frame(component=c("additive","nonadd","maternal"),
    lower=c(quantile(comp[,add],ql_add),quantile(comp[,nonadd],ql_na),quantile(comp[,mater],ql_mat)),
    median= c(quantile(comp[,add],md_add),quantile(comp[,nonadd],md_na),quantile(comp[,mater],md_mat)),
    upper= c(quantile(comp[,add],qu_add),quantile(comp[,nonadd],qu_na),quantile(comp[,mater],qu_mat)))
  ci2<- cbind(ci2$component,round(ci2[,2:4],rnd_r));colnames(ci2)[1]<- "component"
  rownames(ci2)<- 1:3
  ci2_p<- data.frame(component=c("additive","nonadd","maternal"),
    lower=c(quantile(comp$p_add,ql_add),quantile(comp$p_na,ql_na),quantile(comp$p_mat,ql_mat)),
    median= c(quantile(comp$p_add,md_add),quantile(comp$p_na,md_na),quantile(comp$p_mat,md_mat)),
    upper= c(quantile(comp$p_add,qu_add),quantile(comp$p_na,qu_na),quantile(comp$p_mat,qu_mat)))
  ci2_p<- cbind(ci2_p$component,round(ci2_p[,2:4],rnd_p));colnames(ci2_p)[1]<- "component"
  rownames(ci2_p)<- 1:3
if (z0_add == Inf | z0_add == -Inf | z0_na == Inf | z0_na == -Inf | z0_mat == Inf | z0_mat == -Inf)
  { ci2$change<- NA; ci2_p$change<- NA  }
if (z0_add == Inf | z0_add == -Inf) { ci2[1,]<- ci[1,];ci2$change[1] <- "bias fail"
  ci2_p[1,]<- ci_p[1,];ci2_p$change[1] <- "bias fail" }
if (z0_na == Inf | z0_na == -Inf) { ci2[2,]<- ci[2,];ci2$change[2] <- "bias fail"
  ci2_p[2,]<- ci_p[2,];ci2_p$change[2] <- "bias fail" }
if (z0_mat == Inf | z0_mat == -Inf) { ci2[3,]<- ci[3,];ci2$change[3] <- "bias fail"
  ci2_p[3,]<- ci_p[3,];ci2_p$change[3] <- "bias fail" }
} #end ci adjusted
  if (is.null(trait) == T && is.null(bias)) { ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
  if (is.null(trait) == T && !is.null(bias)) { ci_obj<- list(raw=ci2,percentage=ci2_p); return(ci_obj) }
  if (is.null(trait) == F && is.null(bias)) { ci$trait<- as.factor(trait); ci_p$trait<- as.factor(trait)
    ci_obj<- list(raw=ci,percentage=ci_p); return(ci_obj) }
  if (is.null(trait) == F && !is.null(bias)) { ci2$trait<- as.factor(trait); ci2_p$trait<- as.factor(trait)
    ci_obj<- list(raw=ci2,percentage=ci2_p); return(ci_obj) }
}
