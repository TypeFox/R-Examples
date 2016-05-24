blvcm_bin <-
function(pheno, geno, model=3, iter = 30000, burnin= 500, var = -1, lambda = 0.2, cov = 0, init = c(0,0))
{

if ( ( ! is.matrix(pheno) ) || ( ! is.matrix(geno) ) ) {
stop("Phenotypic and genotypic data must be matrices!")
}

num_sub <- nrow(pheno)
num_rv <- ncol(geno)
num_cov <- 0

if(is.matrix(cov))
{
num_cov <- ncol(cov)
}
else
{
if(cov==0)
{cov <- matrix(0)}
else
{
stop("Covariate data must be a matrix!")
}
}

if(min(geno %in% c(0,1,-9))==0)
{
stop("Genotypic data must be 0,1 or -9!")
}

if(lambda <= 0)
{
stop("Lambda must be a positive integer!")
}

burnin <- ceiling(burnin)

if(burnin <= 0)
{
stop("Burnin must be a positive integer!")
}

iter <- ceiling(iter)

if(iter <= 0)
{
stop("The number of MCMC iterations must be a positive integer!")
}


if(burnin > iter)
{
stop("The number of MCMC iterations must be larger than burnin!")
}

if(!model %in% c(1,2,3))
{
stop("Model must be 1 (independent), 2 (AE) or 3 (ACE)")
}


## imputation
geno_temp <- geno
geno_temp[geno_temp==-9] <- 0
miss <- apply(geno, 2, function(x) sum(x==-9))
if(max(miss)==num_sub)
{
stop("Some rare variant has missing data for all individuals.")
}
maf <- apply(geno_temp,2,sum)/(num_sub-miss)

for(i in 1:num_rv)
{
	if(miss[i]>0)
	{
		geno[,i] <- ifelse(geno[,i]==-9, rbinom(1,1,maf[i]), geno[,i])
	}
}

check <- apply(geno,2,max)
if(min(check)==0)
{
stop("There is a column with all 0s in the genotypic data!")
}

num_zyg <- sum(pheno[,3]==1)+sum(pheno[,3]==2)
if(num_zyg != num_sub)
{
stop("Twin type must be 1 or 2!")
}

if(min(pheno[,1] %in% c(0,1))==0)
{
stop("Phenotype must be 1 or 0!")
}

# default var=1
if(var==-1)
{ 
sd_prior <- 1
var_beta <- 1
}
else
{
sd_prior <- sqrt(var)
var_beta <- var
}

## negative init
if(init[1]<0)
{
stop("Initial value of beta must be non-negative.")
}


multResult <- rep(0,13+num_rv+num_cov)

output =.C("CWrapper_blvcm_bin",
product = as.double(multResult),
nRows = as.integer(num_sub),
nCols = as.integer(num_rv),
nCols2 = as.integer(num_cov),
matrix1 = as.integer(t(pheno)),
matrix2 = as.integer(t(geno)),
matrix3 = as.double(t(cov)),
arg_v = as.double(var_beta),
arg_i = as.double(lambda),
arg_n = as.integer(iter),
arg_b = as.integer(burnin),
arg_t = as.integer(model),
arg_i_b = as.double(init[1]),
arg_i_g = as.double(init[2])
)

re_rv <- c()
if(num_rv>0)
{re_rv <- output$product[14:(13+num_rv)]}

re_cov <- c()
if(num_cov>0)
{re_cov <- output$product[(14+num_rv):(13+num_rv+num_cov)]}

# default lambda=0.2
if(lambda==-1)
{error_fun <- 2*(pnorm(0.2, 0,sd_prior)-0.5)}
else
{error_fun <- 2*(pnorm(lambda, 0,sd_prior)-0.5)}

prior_odds_num <- (1-error_fun)*(1-(1/3)^num_rv)
prior_odds_den <- error_fun*(1-(1/3)^num_rv)+(1/3)^num_rv
post_odds_ratio <- output$product[5]/(iter-burnin-output$product[5])
if(post_odds_ratio==Inf)
{
post_odds_ratio <- iter-burnin
warning("The estimated Bayes factor of the main effect is the lower bound due to the limitation of iterations.")
}
bf_main <- post_odds_ratio/(prior_odds_num/prior_odds_den)

inter_num <- output$product[7]
prior_odds_num <- 1-error_fun
prior_odds_den <- error_fun
post_odds_ratio_int <- output$product[6]/(iter-burnin-output$product[6])
if(post_odds_ratio_int==Inf)
{
post_odds_ratio_int <- iter-burnin
warning("The estimated Bayes factor of the interaction effect is the lower bound due to the limitation of iterations.")
}
bf_int <- post_odds_ratio_int/(prior_odds_num/prior_odds_den)

com_var_a <- output$product[9]
com_var_c <- output$product[10]
if(model==2)
{
	com_var_c <- NA
}

if(model==1)
{
	com_var_a <- NA
	com_var_c <- NA
}

return(list(BF_main = bf_main, BF_int = bf_int, post_odds_beta = post_odds_ratio, post_odds_gamma = post_odds_ratio_int,
mean_mu = output$product[1], sd_mu = output$product[11], mean_beta = output$product[2], sd_beta = output$product[12], mean_gamma = output$product[3], 
sd_gamma = output$product[13], zero_alpha=output$product[4],
com_a = com_var_a, com_c = com_var_c, mean_rv = re_rv, mean_cov=re_cov))


}
