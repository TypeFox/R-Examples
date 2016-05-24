c_score=function(mean_sig,var_sig,lambda='default'){
  
  ###########combining site level scores##############
  com_stat=mean_sig
  com_stat$mean_score=-qnorm(mean_sig$mean_p)####computing mean score####
  com_stat=com_stat[order(com_stat$cpg),]
  var_sig=var_sig[order(var_sig$cpg),]
  var_p=var_sig$var_p
  com_stat=cbind(com_stat,var_p)
  com_stat$var_score=-qnorm(com_stat$var_p)####computing var score####3
  com_stat$mean_score[which(com_stat$mean_score<=0)]=0
  com_stat$var_score[which(com_stat$var_score<=0)]=0 
  com_stat$weight=com_stat[,4]+com_stat[,6]
  com_stat=com_stat[which(com_stat$weight!=0),]  ###delete cpg sites whose mean and var scores both smaller than 0#####
  if (lambda=="default"){
    com_stat$ratio=com_stat[,6]/(com_stat[,4]+com_stat[,6])###site level ration####
    tapply(com_stat$ratio, com_stat$gene, function(x){return(mean(x))}) -> mean_lambda
    data.frame(gene=names(mean_lambda), mean_lambda=mean_lambda) -> mean_lambda  #####gene level ratio####3
    com_stat=merge(com_stat,mean_lambda,by="gene")
    lambda_default=unique(com_stat[,c(1,9)])
    lambda_default=mean(lambda_default[,2]) ######computing lambda####
    score=as.numeric(com_stat[,4])*lambda_default+as.numeric(com_stat[,6])*(1-lambda_default)
    com_sig=cbind(com_stat[,c(1,2)],score)  ####site level combined score####
  }  else if (lambda>=0 && lambda<=1){
    score=as.numeric(com_stat[,4])*lambda+as.numeric(com_stat[,6])*(1-lambda)
    com_sig=cbind(com_stat[,c(1,2)],score) 
  } else stop("Please input right lambda\n")
  
  
  
  
  ########computing gene leve score######### 
  tapply(com_sig$score, com_sig$gene, max) -> gene2weight
  data.frame(gene=names(gene2weight), weight=gene2weight) -> gene2weight
  gene2weight = gene2weight[!is.na(gene2weight[,2]), ]
  
  return(gene2weight)
  
}