## Fix by
fixGenotype <- function(genotype){
  for(g in 1:ncol(genotype)){
    afreq <- mean(genotype[,g], na.rm=TRUE)
    if(!is.na(afreq) && afreq > 0.5)
      genotype[,g] <- 2-genotype[,g]
  }
  return(genotype)
}
#genotype <- cbind(rbinom(100,2,0.1), rbinom(100,2,0.9), rep(0,100), rep(1,100), rep(2,100), rep(NA,100))
#print(genotype)
#print(fixGenotype(genotype))

rareGeneTest <- function(genotype, phenotype, use_sign=TRUE, use_weight=TRUE, nperm=1000, binary=all(phenotype==0 | phenotype==1, na.rm=TRUE), strategy="step", thresh=1){
  genotype <- fixGenotype(genotype)

  if(nrow(genotype) != length(phenotype))
    stop("genotype must have the same number of rows as the length of case.")

  p_use_weight = 0
  if(binary & use_weight)
    p_use_weight = 1
  if(!binary & use_weight)
    p_use_weight = 3

  # strategy = 1 (HARD THRESH), 2 (ALL AFREQ), 3 (STEP UP), 4 (ALL SUBSETS)
  p_strategy = 1
  if(strategy == "thresh"){
    p_strategy = 1
  }else if(strategy == "afreq"){
    p_strategy = 2
  }else if(strategy == "step"){
    p_strategy = 3
  }else{
    stop("'strategy' should be one of 'thresh', 'afreq', or 'step'.")
  }

  pvalue = zstat_perm(g=as.matrix(genotype), aff=phenotype, use_sign=use_sign, use_weight=p_use_weight, strategy=p_strategy, nperm=nperm, thresh=thresh)
  return(pvalue)
}

rarePathwayTest <- function(genotype, genotype_gene, phenotype, use_sign=TRUE, use_weight=TRUE, nperm=1000, binary=all(phenotype==0 | phenotype==1, na.rm=TRUE), CUT=15){
  genotype = fixGenotype(genotype)


  t = table(genotype_gene)
  for(ti in 1:length(t)){
    if(t[ti] > CUT){
      wh = which(genotype_gene == names(t)[ti])
      genotype_gene[wh] = paste(genotype_gene[wh], floor(seq(from=1, to=ceiling(t[ti]/15)+0.999999, length.out=t[ti])))
    }
  }


  if(is.character(genotype_gene))
    genotype_gene = as.numeric(as.factor(genotype_gene))

  if(nrow(genotype) != length(phenotype))
    stop("genotype must have the same number of rows as the length of case.")

  p_use_weight = 0
  if(binary & use_weight)
    p_use_weight = 1
  if(!binary & use_weight)
    p_use_weight = 3

  pvalue = zstat_pathway_perm(g=as.matrix(genotype), m=genotype_gene, aff=phenotype, use_sign=use_sign, use_weight=p_use_weight, strategy=3, nperm=nperm)
  return(pvalue)
}
