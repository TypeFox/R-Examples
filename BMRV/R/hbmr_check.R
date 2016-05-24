hbmr_check <- function(geno, qi, fam, kin, iter, burnin, imp, num_sub, num_rv, num_cov)
{

if((sum(geno>1)>0) | (sum(geno<0)>0))
{
stop("Genotypic data must be 0 or 1!")
}

if((imp<=0) | (imp>=1))
{
stop("Imp has to be a real number in (0, 1)!")
}

if(burnin < 0)
{
stop("Burnin must be a positive integer!")
}

if(iter < 0)
{
stop("The number of MCMC iterations must be a positive integer!")
}


if(burnin > iter)
{
stop("The number of MCMC iterations must be larger than burnin!")
}

check <- apply(geno,2,max)
if(min(check)==0)
{
stop("There is a column with all 0s in the genotypic data!")
}

if(!((fam==0) | (fam==1)))
{
stop("The parameter 'fam' must be 0 or 1!")
}

}