# Functions for printing genotypes during the MCMC analysis. Each individual gets their
# own file.

print_g<-function(k,geno.mat,tot.mat,geno_dir){
  rnames<-row.names(tot.mat)
  for(i in 1:nrow(tot.mat)){
    cat(k,geno.mat[i,],file=paste("./",geno_dir,"/",rnames[i],"_g-mcmc.out",sep=""),append=TRUE)
    cat("\n",file=paste("./",geno_dir,"/",rnames[i],"_g-mcmc.out",sep=""),append=TRUE)
  }
}
