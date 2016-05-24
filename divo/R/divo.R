#divo R module 
#Version: 0.1.2
#Autor: Maciej Pietrzak, Michal Seweryn, Grzegorz Rempala
#Maintainer: Maciej Pietrzak <pietrzak.20@osu.edu>
#License: GPL (>=2)

cvg<-function(x){y=x; n = sum(y);f1 = sum(y == 1); if (f1 == n){f1 = n - 1}
  return(1 - f1/n)}

i.inp <- function(x, alpha=1, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if (alpha>0){if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	if (0<CI & CI<1){npySave("tmp_in.pyc", xx, mode="w"); path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/")
  cmd=paste(c("python"), path, c("INP"), alpha,  resample, CI, f='x', PlugIn, format(size, scientific=FALSE), "beta", CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc");    
 	return (.e.clust(xx, csv_output, graph, funct="INP", CVG, PlugIn))}
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}
	else{stop("alpha must be greater than 0")}}}

li <- function(x, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	if (0<CI & CI<1){npySave("tmp_in.pyc", xx, mode="w"); path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/")                                                                                                    
  cmd=paste(c("python"), path, c("LI"), 1,  resample, CI, f='x',PlugIn, format(size, scientific=FALSE), 'beta', "CVG", saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc");  
  return (.e.clust( xx, csv_output, graph, funct="LI", F, PlugIn))}                                                                                     
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}}

ji <- function(x, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, saveBootstrap=FALSE) 
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	if (0<CI & CI<1){npySave("tmp_in.pyc", xx, mode="w"); path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/") 
  cmd=paste(c("python"), path, c("JI"), 1,  resample, CI, f='x',PlugIn, format(size, scientific=FALSE), 'beta', "CVG", saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc");
  return (.e.clust( xx, csv_output, graph, funct="JI", F, PlugIn))}
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}}

rd <- function(x, alpha=0.5, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	chck=0; n=ncol(xx);
	for(a in seq(1,n,1)){if(a<n){for (b in seq(a+1,n,1)){if(sum(x[,a]*x[,b])==0){chck=chck+1;
	warning(as.character(paste("Populations: ", as.character(colnames(x)[a]),", ", as.character(colnames(x)[b]), " are orthogonal.", sep=""))) }}}};
	if(chck==0){if (alpha>0){if (0<CI & CI<1){npySave("tmp_in.pyc", xx, mode="w"); 
	path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/")
  cmd=paste(c("python"), path, c("RD"), alpha,  resample, CI, f='RE',PlugIn, format(size, scientific=FALSE), 'beta', CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc");  
  return (.e.clust( xx, csv_output, graph, funct="RD", CVG, PlugIn))}
	else{stop("Confidence interval must be between 0 and 1; default CI=0.95")}}
	else{stop("alpha must be greater than 0")}}
	else{stop("alpha must be != 1")}}}

srd <- function(x, alpha=0.5, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE)
  {if(alpha < 1){xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if(class(xx)!="matrix"){ xx<-as.matrix(xx)}		
	chck=0; n=ncol(xx);
	for(a in seq(1,n,1)){if(a<n){for (b in seq(a+1,n,1)){if(sum(x[,a]*x[,b])==0){chck=chck+1;
	warning(as.character(paste("Populations: ", as.character(colnames(x)[a]),", ", as.character(colnames(x)[b]), " are orthogonal.", sep=""))) }}}};
	if(chck==0){if (alpha>0){if (0<CI & CI<1){if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	npySave("tmp_in.pyc", xx, mode="w"); path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/")
  cmd=paste(c("python"), path, c("RDS"), alpha,  resample, CI, f='RE',PlugIn, format(size, scientific=FALSE), 'beta', CVG, saveBootstrap, sep=" "); t<-system(cmd, intern=TRUE); 
  file.remove("tmp_in.pyc");  
  return (.e.clust( xx, csv_output, graph, funct="RD", CVG, PlugIn))}
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}
	else{stop("alpha must be greater than 0")}}
	else {stop("alpha must be < 1")}}}}

mh <- function(x, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if(0<CI & CI<1){
	if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	npySave("tmp_in.pyc", xx, mode="w"); path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/") 
  cmd=paste(c("python"), path, c("MH"), 1,  resample, CI, f='RE',PlugIn, format(size, scientific=FALSE), 'beta', "CVG", saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc"); 
  return (.e.clust( xx, csv_output, graph, funct="MH", F, PlugIn))}       
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}}

pg <- function(x, alpha=1, beta=alpha, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if (alpha>0 & beta>0){if (0<CI & CI<1){if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	npySave("tmp_in.pyc", xx, mode="w"); path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/")
  cmd=paste(c("python"), path, c("PG"), alpha,  resample, CI, f='x',PlugIn, format(size, scientific=FALSE), beta, CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc"); 
  return (.e.clust( xx, csv_output, graph, funct="PG", CVG, PlugIn))}
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}
	else{stop("alpha and beta must be greater than 0")}}}

pg.ht <- function(x, alpha=1,  beta=alpha, CI=0.95, resample=100, graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if (alpha>0 & beta>0){if (0<CI & CI<1){if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	npySave("tmp_in.pyc", xx, mode="w"); path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/")
  cmd=paste(c("python"), path, c("PG_HT"), alpha,  resample, CI, f='x',PlugIn, format(size, scientific=FALSE), beta, CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc"); 
  return (.e.clust( xx, csv_output, graph, funct="PG.ht", CVG, PlugIn))}
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}
	else{stop("alpha and beta must be greater than 0")}}}

i.in <- function(x, alpha=1, CI=0.95, resample=100, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  	wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
	if(class(xx)=="numeric"){stop("Number of columns must be greater than 1")}
	if(resample < 1){stop("resample must be greater than 1")}
	if(size < 0.1){stop("size must be greater than or equal to 0.1")}
	if (0.1<=alpha & alpha<=1) {xx=x; if (0<CI & CI<1){
	if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
	npySave("tmp_in.pyc", xx, mode="w");                                               
  path<-paste(system.file(package="divo"), "DivO_Overlap.py", sep="/")
  cmd=paste(c("python"), path, c("IN"), alpha,  resample, CI, f='x',PlugIn, format(size, scientific=FALSE), 1, CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc"); get.data<-t(npyLoad('tmp_IN.npy'));
  rownames(get.data)<-c("Mean", "Lower.Quantile", "Upper.Quantile"); file.remove("tmp_IN.npy");
  return (get.data)}
	else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}}
	else{stop("alpha must be between 0.1 and 1")}}}

.e.clust<-function( xy, csv_output, graph, funct, CVG, PlugIn)
  {pr_mean<-t(npyLoad('tmp_out_mean.npy')); pr_min <-t(npyLoad('tmp_out_min.npy')); pr_max <-t(npyLoad('tmp_out_max.npy'));
  colnames(pr_mean) <- colnames(xy) -> rownames(pr_mean); colnames(pr_min) <- colnames(xy) -> rownames(pr_min); colnames(pr_max) <- colnames(xy) -> rownames(pr_max);
  try(if(funct !="RD"){diag(pr_mean)<-1;diag(pr_min)<-1;diag(pr_max)<-1}, silent=TRUE)
  outlist<-list(pr_mean, pr_min, pr_max); names(outlist) <- c("Mean", "Lower.Quantile", "Upper.Quantile") 
  if(graph!=TRUE & graph!=FALSE){fname=paste(graph, '.pdf', sep='')}; if(graph==TRUE){fname="DivO_Clust.pdf"}
  if(graph!=FALSE){.p.clust(xy, fname, funct, lab=xy[0,], PlugIn)}  
  if(csv_output!=FALSE & csv_output!=TRUE){fname.csv=paste(csv_output, '.csv',sep='')}; if(csv_output==TRUE){fname.csv="DivO_Overlap_out.csv"}
  if(csv_output!=FALSE){f.list.to.df<-function(csv.x){colnames(csv.x[[1]])->cnames; colnames(csv.x[[1]])<-NULL; rownames(csv.x[[1]])<-NULL
  temp<-rep(0,length(csv.x[[1]][1,]));for(i in 1:length(csv.x)){colnames(csv.x[[i]])<-NULL;rownames(csv.x[[i]])<-NULL;temp<-rbind(temp,rep('',length(csv.x[[1]][1,])),csv.x[[i]])};temp<-temp[-1,]
  temp<-cbind(c('Mean',cnames,'Lower.Quantile',cnames,'Upper.Quantile',cnames),temp); colnames(temp)<-c('Type',cnames)
  temp[temp==0]<-''; data.frame(temp)->temp; rownames(temp)<-NULL
  return(temp)}  
  f.list.to.df(outlist)->overlap.list.to.df
  write.csv(file=fname.csv, x=overlap.list.to.df, row.names=FALSE)};
  file.remove("tmp_out_mean.npy", "tmp_out_min.npy", "tmp_out_max.npy")
  if(PlugIn==TRUE){outlist<-outlist[1]; names(outlist)[length(outlist)]<- "PlugIn"}
  if(CVG==TRUE){a=npyLoad("tmp_out_CVG.npy"); a=t(rbind(a)); 
  	colnames(a)<-c("CVG");rownames(a)<-colnames(xy);if(PlugIn==TRUE){colnames(a)<-c("PlugIn.CVG")}
	outlist[[length(outlist)+1]]<-a; names(outlist)[length(outlist)]<- "Coverage"; file.remove("tmp_out_CVG.npy")}
  return (outlist)}

.p.clust<- function(xx,fname, funct, lab, PlugIn, h=FALSE, height=7, width=0.6*height)
  {if(funct != "RD"){
  mat.1<-rbind(lab, 1-as.matrix(round(t(npyLoad("tmp_out_mean.npy")),6))); 
  mat.2<-rbind(lab, 1-as.matrix(round(t(npyLoad("tmp_out_min.npy")),6)));  
  mat.3<-rbind(lab, 1-as.matrix(round(t(npyLoad("tmp_out_max.npy")),6)))}
  else{
  mat.1<-rbind(lab, as.matrix(round(t(npyLoad("tmp_out_mean.npy")),6))); 
  mat.2<-rbind(lab, as.matrix(round(t(npyLoad("tmp_out_min.npy")),6)));  
  mat.3<-rbind(lab, as.matrix(round(t(npyLoad("tmp_out_max.npy")),6)))}
  
  as.dendrogram(agnes(mat.1,diss=TRUE,method='ward'))->dend.mat.1; as.dendrogram(agnes(mat.2,diss=TRUE,method='ward'))->dend.mat.2; as.dendrogram(agnes(mat.3,diss=TRUE,method='ward'))->dend.mat.3; 
  cex.val=0.4; if(ncol(xx) < 5){cex.val=0.5}; if(ncol(xx) > 30){cex.val=0.3};
  if(PlugIn==FALSE){par(mfrow=c(3,1), cex=cex.val, lwd=0.5, mar=c(10, 4, 2, 0.5)); 
  plot(dend.mat.1,horiz = h,main="Mean"); plot(dend.mat.2,horiz = h,main="Lower Quantile"); plot(dend.mat.3,horiz = h,main="Upper Quantile");
  dev.copy(pdf, fname, height=height, width=width); dev.off()}
  if(PlugIn==TRUE){par(mfrow=c(1,1), mar=c(10, 4, 2, 0.5)); plot(dend.mat.1,horiz = h,main="PlugIn");
  dev.copy(pdf, fname); dev.off()}}

dp <- function(x, alpha=seq(0.1, 2, 0.1), CI=0.95, resample=100, single_graph=FALSE, pooled_graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE) 
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
  if(resample < 1){stop("resample must be greater than 1")}
  if(size < 0.1){stop("size must be greater than or equal to 0.1")}
  if (0<CI & CI<1){f="RE";if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
  npySave("tmp_in.pyc", xx, mode="w"); Alpha.Profile=alpha;
  if(min(Alpha.Profile)<=0){stop("alpha must be larger than 0")}
  if(class(Alpha.Profile)!="character"){npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  else{Alpha.Profile=seq(0.1, 2, 0.1); npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  path<-paste(system.file(package="divo"), "DivO_DP.py", sep="/")  
  cmd=paste(c("python"), path, c("DP"), "alpha",  resample, CI, f, PlugIn, format(size, scientific=FALSE), "Alpha.Profile", CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); 
  file.remove("tmp_in.pyc");if(class(Alpha.Profile)!="character"){file.remove("tmp_alpha.pyc")}
  tmp<-.DP.single.profile(c(colnames(xx)), single_graph, "Diversity Profile:", Alpha.Profile, CVG, PlugIn)
  output<-tmp[[1]]; min.y<-tmp[[2]]; max.y<-tmp[[3]]; names(output)<-colnames(xx);out1<-output
  if(CVG==TRUE){out1<-output; output[[length(output)+1]]<-tmp[[4]]; names(output)[length(output)]<- "Coverage"; file.remove("DP_tmp_CVG.npy") }
  if(length(alpha)==1 & class(Alpha.Profile)!="character"){warning("For graphics alpha must be a vector of lenght greater than 1")}
  if(pooled_graph!=FALSE & length(alpha)>1){.DP.pooled.profile(colnames(xx),  min.y, max.y, "DP ", pooled_graph, Alpha.Profile)}
  if(csv_output!=FALSE){if(csv_output==TRUE){write.csv(file="DivO_DP_output.csv", x= out1, row.names=FALSE)}
  else{write.csv(file=paste(csv_output, ".csv", sep=""), x= out1, row.names=FALSE)}}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output')); file.remove(list_to_remove)
  return(output)}     
  else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output_'))
  file.remove(list_to_remove)}}

ens <- function(x, alpha=seq(0.1, 2, 0.1), CI=0.95, resample=100, single_graph=FALSE, pooled_graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE) 
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
  if(resample < 1){stop("resample must be greater than 1")}  	
  if(size < 0.1){stop("size must be greater than or equal to 0.1")}
  if (0<CI & CI<1){
  f="RE"
  if(class(xx)!="matrix"){ xx<-as.matrix(xx)} 
  npySave("tmp_in.pyc", xx, mode="w"); 
  Alpha.Profile=alpha
  if(min(Alpha.Profile)<=0){stop("alpha must be larger than 0")}
  if(class(Alpha.Profile)!="character"){npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  else{Alpha.Profile=seq(0.1, 2, 0.1); npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  path<-paste(system.file(package="divo"), "DivO_DP.py", sep="/")
  cmd=paste(c("python"), path, c("ENS"), "alpha",  resample, CI, f, PlugIn, format(size, scientific=FALSE), "Alpha.Profile", CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); 
  file.remove("tmp_in.pyc"); 
  if(class(Alpha.Profile)!="character"){file.remove("tmp_alpha.pyc")} 
  tmp<-.DP.single.profile(c(colnames(xx)), single_graph, "Diversity Profile:", Alpha.Profile, CVG, PlugIn)
  output<-tmp[[1]]; min.y<-tmp[[2]]; max.y<-tmp[[3]];names(output)<-colnames(xx);out1<-output
  if(CVG==TRUE){out1<-output; output[[length(output)+1]]<-tmp[[4]]; names(output)[length(output)]<- "Coverage"; file.remove("DP_tmp_CVG.npy") }
  if(length(alpha)==1 & class(Alpha.Profile)!="character"){warning("For graphics alpha must be a vector of lenght greater than 1")}
  if(pooled_graph!=FALSE & length(alpha)>1){.DP.pooled.profile(colnames(xx), min.y, max.y, "ENS", pooled_graph, Alpha.Profile)}
  if(csv_output!=FALSE){if(csv_output==TRUE){write.csv(file="DivO_DP_output.csv", x= out1, row.names=FALSE)}
  else{write.csv(file=paste(csv_output, ".csv", sep=""), x= out1, row.names=FALSE)}}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output'))
  file.remove(list_to_remove)
  return(output)}     
  else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output_'))
  file.remove(list_to_remove)}}

dp.ht <- function(x, alpha=seq(0.1, 2, 0.1), CI=0.95, resample=100, single_graph=FALSE, pooled_graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE) 
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
  if(resample < 1){stop("resample must be greater than 1")}
  if(size < 0.1){stop("size must be greater than or equal to 0.1")}
  if (0<CI & CI<1){
  if(class(xx)!="matrix"){ xx<-as.matrix(xx)}
  Alpha.Profile=alpha
  if(min(Alpha.Profile)<=0){stop("alpha must be larger than 0")}
  npySave("tmp_in.pyc", xx, mode="w"); f="HT";
  if(class(Alpha.Profile)!="character"){npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  else{Alpha.Profile=seq(0.1, 2, 0.1); npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  path<-paste(system.file(package="divo"), "DivO_DP.py", sep="/")
  cmd=paste(c("python"), path, c("DP"), "alpha",  resample, CI, f, PlugIn, format(size, scientific=FALSE), "Alpha.Profile", CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); 
  file.remove("tmp_in.pyc"); 
  if(class(Alpha.Profile)!="character"){file.remove("tmp_alpha.pyc")}
  tmp<-.DP.single.profile(c(colnames(xx)), single_graph, "Diversity Profile:", Alpha.Profile, CVG, PlugIn)
  output<-tmp[[1]]; min.y<-tmp[[2]]; max.y<-tmp[[3]];names(output)<-colnames(xx);out1<-output
  if(CVG==TRUE){out1<-output; output[[length(output)+1]]<-tmp[[4]]; names(output)[length(output)]<- "Coverage"; file.remove("DP_tmp_CVG.npy") }
  if(length(alpha)==1 & class(Alpha.Profile)!="character"){warning("For graphics alpha must be a vector of lenght greater than 1")}
  if(pooled_graph!=FALSE & length(alpha)>1){.DP.pooled.profile(colnames(xx),  min.y, max.y, "DP ", pooled_graph, Alpha.Profile)}                   
  if(csv_output!=FALSE){if(csv_output==TRUE){write.csv(file="DivO_DP_output.csv", x= out1, row.names=FALSE)}
  else{write.csv(file=paste(csv_output, ".csv", sep=""), x= out1, row.names=FALSE)}}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output'))
  file.remove(list_to_remove)
  return(output)}     
  else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output_'))
  file.remove(list_to_remove)}}

ens.ht <- function(x, alpha=seq(0.1, 2, 0.1), CI=0.95, resample=100, single_graph=FALSE, pooled_graph=FALSE, csv_output=FALSE, PlugIn=FALSE, size=1, CVG=FALSE, saveBootstrap=FALSE)
  {xx=x; resample= round(resample, 0); xx[is.na(xx)] <- 0;
  wdCheck<-try(.WriteTest(), silent=T);if(wdCheck=="Passed"){
  if(resample < 1){stop("resample must be greater than 1")}
  if(size < 0.1){stop("size must be greater than or equal to 0.1")}
  if (0<CI & CI<1){
  f="HT";
  if(class(xx)!="matrix"){ xx<-as.matrix(xx)} 
  npySave("tmp_in.pyc", xx, mode="w"); 
  Alpha.Profile=alpha
  if(min(Alpha.Profile)<=0){stop("alpha must be larger than 0")}
  if(class(Alpha.Profile)!="character"){npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  else{Alpha.Profile=seq(0.1, 2, 0.1); npySave("tmp_alpha.pyc", Alpha.Profile, mode="w")}
  path<-paste(system.file(package="divo"), "DivO_DP.py", sep="/")
  cmd=paste(c("python"), path, c("ENS"), "alpha",  resample, CI, f, PlugIn, format(size, scientific=FALSE), "Alpha.Profile", CVG, saveBootstrap, sep=" "); 
  t<-system(cmd, intern=TRUE); file.remove("tmp_in.pyc"); 
  if(class(Alpha.Profile)!="character"){file.remove("tmp_alpha.pyc")} 
  tmp<-.DP.single.profile(c(colnames(xx)), single_graph, "Diversity Profile:", Alpha.Profile, CVG, PlugIn)
  output<-tmp[[1]]; min.y<-tmp[[2]]; max.y<-tmp[[3]];names(output)<-colnames(xx);out1<-output
  if(CVG==TRUE){out1<-output; output[[length(output)+1]]<-tmp[[4]]; names(output)[length(output)]<- "Coverage"; file.remove("DP_tmp_CVG.npy") }
  if(length(alpha)==1 & class(Alpha.Profile)!="character"){warning("For graphics alpha must be a vector of lenght greater than 1")}
  if(pooled_graph!=FALSE & length(alpha)>1){.DP.pooled.profile(colnames(xx), min.y, max.y, "ENS", pooled_graph, Alpha.Profile)}
  if(csv_output!=FALSE){if(csv_output==TRUE){write.csv(file="DivO_DP_output.csv", x= out1, row.names=FALSE)}
  else{write.csv(file=paste(csv_output, ".csv", sep=""), x= out1, row.names=FALSE)}}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output')); file.remove(list_to_remove)
  return(output)}     
  else {stop("Confidence interval must be between 0 and 1; default CI=0.95")}
  list_to_remove<-(list.files(pattern= 'DP_tmp_output_'))
  file.remove(list_to_remove)}}

.DP.single.profile<-function(lab, single_graph, title='Title', Alpha.Profile, CVG, PlugIn)
  {max.y<-0; min.y<-100; pat<-c('DP_tmp_output_');  allfiles<- list.files(pattern= pat); 
  
  ###
  #len<-seq(1,length(allfiles), by=1); 
  len<-seq_along(allfiles)
  ###
  if(length(lab)<=1){len=1}
  xx_pooled<-list() 
  if(CVG==TRUE){cvg.tmp<-npyLoad("DP_tmp_CVG.npy");cvg.tmp=cbind(cvg.tmp);rownames(cvg.tmp)<-lab; colnames(cvg.tmp)<-"CVG"; if(PlugIn==TRUE){colnames(cvg.tmp)<-"PlugIn.CVG"}}
  for (l in len){ 
  fname<-allfiles[l]
  xx<-npyLoad(fname, dotranspose=TRUE);
  if(sum(xx[,1]-Alpha.Profile)!=0){xx<-npyLoad(fname, dotranspose=FALSE)}
  mmat<-xx[,2]; q1<-xx[,2]-xx[,3];
  q2<-xx[,2]+xx[,4];alfa<-xx[,1];
  sv<-cbind(xx[,1],xx[,2],q1,q2);
  if(max(q2) > max.y){max.y<-max(q2)}; 
  if(min(q1) < min.y){min.y<-min(q1)};
  y.label='Diversity Index'
  if(title=="ENS:"){y.label="Effective Number of Species"}
  if(single_graph!=FALSE){par(mfrow=c(1,1))
  plot(alfa,mmat,type='l',ylim=c(min(q1),max(q2)), main= paste(title, sub(".*put_", "",  sub(".npy", "", lab[l]))), xlab='Order (alpha)',ylab=y.label, cex.axis=0.8, cex.lab=0.9, cex.main=0.9)  
  lines(alfa,mmat,type='l', col='black'); lines(alfa,q1,lty=3, col='black'); lines(alfa,q2,lty=3, col='black') 
  if(single_graph==TRUE){dev.copy(pdf, paste("DivO_SingleProfile_", sub(".*put_", "", substr(allfiles[l], 1, nchar(allfiles[l])-4)), '.pdf', sep="")); dev.off() }
  else{dev.copy(pdf, paste(single_graph, "_",sub(".*put_", "", substr(allfiles[l], 1, nchar(allfiles[l])-4)), '.pdf', sep="")); dev.off() }}
  colnames(sv)<-c('Alpha', 'Mean', 'Lower.Quantile', 'Upper.Quantile')
  if(PlugIn==TRUE){sv=sv[,1:2]; colnames(sv)<-c("Alpha", "PlugIn")}
  xx_pooled[[length(xx_pooled)+1]]<-sv}
  if(CVG==TRUE){l=list(xx_pooled, min.y, max.y, cvg.tmp)}
  else{l=list(xx_pooled, min.y, max.y)}
  return(l)}

.DP.pooled.profile<-function( lab, min.y, max.y, title="Title", pooled_graph, Alpha.Profile)
  {pat<-c('DP_tmp_output_'); allfiles<- list.files(pattern= pat); 
  	###
  	###len<-seq(1,length(allfiles), by=1); 
  	len<-seq_along(allfiles)
  	###
  if(length(lab)==1){len=1}
  xx_pooled<-list(); colors<-rainbow(length(lab));
  par(mfrow=c(1,1)); y.label='Diversity Index'
  if(title=="ENS"){y.label="Effective Number of Species"}
  fname=allfiles[1]; xx<-npyLoad(fname, dotranspose=TRUE);
  if(sum(xx[,1]-Alpha.Profile)!=0){xx<-npyLoad(fname, dotranspose=FALSE)}
  mmat<-xx[,2]; q1<-xx[,2]-xx[,3] ;q2<-xx[,2]+xx[,3];alfa<-xx[,1];
  plot(alfa,mmat,type='l',ylim=c(min.y, max.y), main=paste(title, 'All Populations'), xlab='Order (alpha)',ylab=y.label, cex.axis=0.8, cex.lab=0.9, cex.main=0.9)      
  for (l in len){fname=allfiles[l];
  xx<-npyLoad(fname, dotranspose=TRUE);
  if(sum(xx[,1]-Alpha.Profile)!=0){xx<-npyLoad(fname, dotranspose=FALSE)}
  mmat<-xx[,2]; q1<-xx[,2]-xx[,3] ; q2<-xx[,2]+xx[,3]; alfa<-xx[,1]; 
  lines(alfa,mmat,type='l', col=colors[l]); lines(alfa,q1,lty=3, col=colors[l]); lines(alfa,q2,lty=3, col=colors[l])}  
  list.xx<-lab; 
  try(legend('topright', list.xx, cex=0.5, col=c(colors[1:length(list.xx)]), lty=1:1, lwd=1, bty = "n", xpd=NA), silent=TRUE);
  if(pooled_graph==TRUE){ dev.copy(pdf, paste("DivO_DP_AllPopulations", '.pdf', sep="")); dev.off()}
  else{dev.copy(pdf, paste(pooled_graph, '.pdf', sep="")); dev.off()}
  return()}
  
.WriteTest <- function(){testMatrix<-NULL->np;
  path<-paste(system.file(package="divo"), "DivO_EnvCheck.py", sep="/");
  try(npySave("WriteTest", matrix(1,2,2)), silent=T)
  try(testMatrix<-npyLoad("WriteTest"), silent=T)
  cmd=paste(c("python"), path, sep=" "); 
  try(t<-system(cmd, intern=TRUE), silent=T); 
  try(np<-npyLoad("divoEnvCheck_out.npy"), silent=T);
  options(warn=-1);try(file.remove("WriteTest", "divoEnvCheck_out.npy"), silent=T)
  if(sum(testMatrix)==4 & sum(np)==12){return("Passed")}
  else {message("\ndivo environment check: Failed. \n\nSolutions:\n1. check write permissions for working directory,\n2. check if dependences are installed and configured correctly.\n");return("Failed")}}

