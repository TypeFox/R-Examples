permutation <- function (normal_data,tumor_data,annotation,res,n=100){
  
  c_score= function (mean_sig, var_sig, lambda = "default") 
  {
    com_stat = mean_sig
    com_stat$mean_score = -qnorm(mean_sig$mean_p)
    com_stat = com_stat[order(com_stat$cpg), ]
    var_sig = var_sig[order(var_sig$cpg), ]
    var_p = var_sig$var_p
    com_stat = cbind(com_stat, var_p)
    com_stat$var_score = -qnorm(com_stat$var_p)
    com_stat$mean_score[which(com_stat$mean_score <= 0)] = 0
    com_stat$var_score[which(com_stat$var_score <= 0)] = 0
    com_stat$weight = com_stat[, 4] + com_stat[, 6]
    com_stat = com_stat[which(com_stat$weight != 0), ]
    if (lambda == "default") {
      com_stat$ratio = com_stat[, 6]/(com_stat[, 4] + com_stat[, 
                                                               6])
      mean_lambda <- tapply(com_stat$ratio, com_stat$gene, 
                            function(x) {
                              return(mean(x))
                            })
      mean_lambda <- data.frame(gene = names(mean_lambda), 
                                mean_lambda = mean_lambda)
      com_stat = merge(com_stat, mean_lambda, by = "gene")
      lambda_default = unique(com_stat[, c(1, 9)])
      lambda_default = mean(lambda_default[, 2])
      score = as.numeric(com_stat[, 4]) * lambda_default + 
        as.numeric(com_stat[, 6]) * (1 - lambda_default)
      com_sig = cbind(com_stat[, c(1, 2)], score)
    }
    else if (lambda >= 0 && lambda <= 1) {
      score = as.numeric(com_stat[, 4]) * lambda + as.numeric(com_stat[, 
                                                                       6]) * (1 - lambda)
      com_sig = cbind(com_stat[, c(1, 2)], score)
    }
    else stop("Please input right lambda\n")
    gene2weight <- tapply(com_sig$score, com_sig$gene, max)
    gene2weight <- data.frame(gene = names(gene2weight), weight = gene2weight)
    gene2weight = gene2weight[!is.na(gene2weight[, 2]), ]
    return(gene2weight)
  }
  
  
  for (j in 1:n){

#############   permuting sick label############3 
    pair_data=data.frame(cbind(normal_data[,1],tumor_data[,1]))
    names(pair_data)=sample(c("N","T"))
    tumor_perm=data.frame(pair_data$T)
    normal_perm=data.frame(pair_data$N)
    for (i in 2:length(normal_data[1,])){
      pair_data=data.frame(cbind(cbind(normal_data[,i],tumor_data[,i])))
      names(pair_data)=sample(c("N","T"))
      tumor_perm=cbind(tumor_perm,pair_data$T)
      normal_perm=cbind(normal_perm,pair_data$T)
    }
  
  
###########site p values########
  meanp_perm_r=rep(999,times=length(normal_data[,1]))
  varp_perm_r=rep(999,times=length(normal_data[,1]))
  for (k in 1:length(normal_data[,1])){
    ttest=t.test(as.numeric(normal_perm[k,1:length(normal_data[1,])]),as.numeric(tumor_perm[k,1:length(normal_data[1,])]),paired=T)
    if (ttest$statistic<0){meanp_perm_r[k]=pt(ttest$statistic, ttest$parameter)*2}  else 
    {meanp_perm_r[k]= 2 * (pt(-ttest$statistic, ttest$parameter))}
    vartest=var.test(as.numeric(normal_perm[k,]),as.numeric(tumor_perm[k,]),paired=T)
    if (vartest$statistic<0){varp_perm_r[k]=pt(vartest$statistic, vartest$parameter)}  else 
    {varp_perm_r[k]= 1}
  }

  
meansig_perm=cbind(annotation,meanp_perm_r)
varsig_perm=cbind(annotation,varp_perm_r)
names(meansig_perm)=c("cpg","gene","mean_p")
names(varsig_perm)=c("cpg","gene","var_p")


weight_perm=c_score(meansig_perm,varsig_perm)
save(weight_perm,file=paste("weight_perm_",j,".RData",sep = ""))
}

perm_score=weight_perm
names(perm_score)[2]=0
for (i in 1:n){
  load(paste("weight_perm_",i,".RData",sep = ""))
  names(weight_perm)[2]=i
  perm_score=merge(perm_score,weight_perm,by="gene")
}
perm_score=perm_score[,-2]

len = floor(length(res$zi.order[, 1]) * 0.01)
idx1 = as.character(names(res$genesets.clear))
idx2 = as.character(res$zi.ordered$gene[1:(len)])
idx = match(idx2, idx1)
module.list= res$genesets.clear[idx]
perm_result=res$zi.ordered[1:(len),1:2]
for (j in 2:n+1){
  perm=rep(999,len)
  for (i in 1:len){
    genes = module.list[[i]]
    match(genes, perm_score[,1]) -> idx3
    idx3 = idx3[!is.na(idx3)]
    perm[i]=sum(as.numeric(perm_score[idx3, j]))/sqrt(length(idx3))
  }
  perm_result=data.frame(cbind(perm_result,perm))
}
for (i in 1:n){
  perm_result$p[i]=length(perm_result[i,3:n+2][which(perm_result[i,3:n+2]>perm_result[i,2])])/100
}

return(perm_result)

}



