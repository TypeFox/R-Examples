#berend krebs master thesis
#math background module



fstcal <- function(i,j){
	return (exp(BBB2$d_alpha[i]+BBB2$d_beta[j])/(1+exp(BBB2$d_alpha[i]+BBB2$d_beta[j])))
}


## log_prior_alpha # Returns Log(pi(alpha))
log_prior_alpha <- function(alpha_){

#return(log(dsn(alpha_, location=-2, scale=0.5, shape=5)))

	return(log(0.5*(1/(BBB2$sd_prior_alpha*sqrt(2*pi)))*exp(-(alpha_-BBB2$m1_prior_alpha)*(alpha_-BBB2$m1_prior_alpha)/(2*BBB2$sd_prior_alpha^2)) + 0.5*(1/   (BBB2$sd_prior_alpha*sqrt(2*pi)))*exp(-(alpha_-BBB2$m2_prior_alpha)*(alpha_-BBB2$m2_prior_alpha)/(2*BBB2$sd_prior_alpha^2))))

}

 

 allelecount_loglikelihood <- function(x){

  x         <- x[BBB2$GROUP]
  theta     <- outer(x,BBB2$d_beta,"+")
  theta     <- exp(-(theta))
  L1        <- lfactorial(BBB2$sample_size)  + lgamma(theta)-lgamma(BBB2$sample_size + theta)
  
  val <- matrix(0,length(BBB2$d_alpha),dim(BBB2$freq_locus)[2])
  for(xx in 1:BBB2$popnum){
    xyz <- theta[,xx]*BBB2$freq_locus
    val <- val + lgamma(BBB2$freq_pop[[xx]] + xyz) - lfactorial(BBB2$freq_pop[[xx]]) - lgamma(xyz)
  }              
  
  return(sum(L1) + sum(val))

   # Original C++
    # theta=exp(-(alpha[i]+beta[j]));
     #               loglikelihood+=factln(pop[j].locus[i].alleleCount)+gammaln(theta)
      #                             -gammaln(pop[j].locus[i].alleleCount+theta);
       #             for (int k=0;k<pop[j].locus[i].ar;k++)
        #                loglikelihood+=gammaln(pop[j].locus[i].data_allele_count[k]+theta*freq_locus[i].allele[k])
         #                              -factln(pop[j].locus[i].data_allele_count[k])-gammaln(theta*freq_locus[i].allele[k]);

}


