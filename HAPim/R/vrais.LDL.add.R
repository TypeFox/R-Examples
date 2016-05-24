`vrais.LDL.add` <-
function(moyenne.pere,alpha.Q,s,CD,perf,PLA,DL.m,DL.chrom1,DL.chrom2,desc.pere,mean.gene)
  {
nb.pere  = length(desc.pere[,1])
logvrais = 0

 for (i in 1:nb.pere)
 {
 deb=desc.pere[i,1]
 fin=desc.pere[i,2]

 vrais.intra.pere = vrais.LDL.add.pere(moyenne.pere[i],alpha.Q,s,CD[deb:fin],perf[deb:fin],PLA[deb:fin],DL.m[deb:fin],DL.chrom1[i],DL.chrom2[i],mean.gene)
 
 logvrais = logvrais + vrais.intra.pere
  }
logvrais

}

