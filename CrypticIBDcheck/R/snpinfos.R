
###############################################
snpinfos=function(snpmatlist){

  csct=snpmatlist$subject.support$csct
  snpobj0=snpmatlist$snp.data[csct==0,] #extracting controls
  Position=snpmatlist$snp.support$Position
  Chromosome=snpmatlist$snp.support$Chromosome
  Gen_loc=snpmatlist$snp.support$Gen_loc


  if(anyDuplicated(cbind(Position,Chromosome))) {
    stop("duplicate physical map positions found in input data set")
  }
    
  max=length(snpobj0[1,])
  sumsnp0<-chopsticks::summary(snpobj0)
  vecpval=2*pnorm(abs(sumsnp0$z.HWE), lower.tail = F)
  Chromosome=gsub("[a-z]", "", Chromosome,perl=TRUE)
   
  if (is.null(Gen_loc)){
    Gen_loc <- SNPgenmap(Position,Chromosome)
  }
  SNP.support=data.frame(Chromosome=as.numeric(Chromosome),
    Position=Position,
    Gen_loc=Gen_loc,
    pvalue_HWE=as.vector(vecpval))
  rownames(SNP.support)=colnames(snpobj0)
  snpmatlist$snp.support=SNP.support


  return(snpmatlist)
}

