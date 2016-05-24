hbmr_ord <-
function(pheno, geno, qi = matrix(), fam = 0, kin = matrix(), iter = 10000, burnin = 500, gq = 20, imp = 0.1, cov = matrix(), maf = c(), pa = 1.3, pb = 0.04)
{

if ( ( ! is.matrix(geno) ) | ( ! is.matrix(qi) ) | ( ! is.matrix(kin) ) | ( ! is.matrix(cov) ) ) {
stop("Genotype, quality, relatedness and covariate data must be matrices!")
}

num_sub <- length(pheno)
num_rv <- ncol(geno)
num_cov <- 0

burnin <- ceiling(burnin)
iter <- ceiling(iter)

hbmr_check(geno, qi, fam, kin, iter, burnin, imp, num_sub, num_rv, num_cov)

if(sum(dim(cov))>2)
{
if(nrow(cov)!=num_sub)
{
stop("Covariate matrix has incorrect rows!")	
}
num_cov <- ncol(cov)
}
else
{
	cov <- matrix(0)
}


G <- matrix(0)
if(fam == 1)
{
if(sum(dim(kin)==c(num_sub,num_sub))<2)
{
stop("The dimension of the kinship matrix is incorrect!")
}
ei_ran <- eigen(kin)
G <- ei_ran$vectors%*%diag(sqrt(ifelse(ei_ran$value<0,0,ei_ran$value)))
}

if((sum(dim(qi))==2) & (is.na(qi[1,1])))
{
	qi <- matrix(99, num_sub, num_rv) 
}

if(length(maf)==0)
{
	maf <- -1
}
else
{
	if(length(maf)!=num_rv)
	{
		stop("The MAF information is incomplete!")
	}
}


y <- pheno
x <- t(geno)
q <- t(qi)
c <- t(cov)
# ran <- matrix(0)
ran <- t(G)

rvinfo <- FALSE 

rv_detail <- 0
if(rvinfo == TRUE)
{rv_detail <- 1}
if(rv_detail == 0)
{
multResult <- rep(0,5+num_cov)
}
else
{
multResult <- rep(0,5+num_cov+2*num_rv)
}

output =.C("CWrapper_hbmr_ord",
product = as.double(multResult),
nRows = as.integer(num_sub),
nCols = as.integer(num_rv),
nCols2 = as.integer(num_cov),
fam = as.integer(fam),
matrix1 = as.integer(y),
matrix2 = as.double(x),
matrix3 = as.double(q),
matrix4 = as.double(c),
kin_m = as.double(ran),
arg_m = as.double(maf),
arg_i = as.double(imp),
arg_n = as.integer(iter),
arg_bu = as.integer(burnin),
arg_t = as.double(gq),
gamma_est = as.integer(rv_detail),
arg_a = as.double(pa),
arg_b = as.double(pb)
)


#return(output$product)

if(output$product[1]!=1)
{
bf <- output$product[1]/(1-output$product[1])
}
else
{
bf <- Inf
}

if(output$product[2]!=1)
{
bf_rb <- output$product[2]/(1-output$product[2])
}
else
{
bf_rb <- iter-burnin
warning("The estimated Bayes factor bf_rb is the lower bound due to the limitation of numerical precision.")
}

re_cov <- c()
if(num_cov>0)
{re_cov <- output$product[6:(5+num_cov)]}

re_rv_es <- c()
if(rvinfo==TRUE)
{
re_rv_es <- output$product[(5+num_cov+1):(5+num_cov+num_rv)]
}

re_rv_sd <- c()
if(rvinfo==TRUE)
{
re_rv_sd <- output$product[(5+num_cov+num_rv+1):(5+num_cov+num_rv+num_rv)]
}

p_val <- NA
if(bf_rb>2)
{
	fun <- function (x) (-1)/(exp(1)*x*log(x))-bf_rb
	p_val <- uniroot(fun, tol = 0.00000000000000000001, interval=c(0.0000000000000000001,0.5))
}

return(list(BF = bf, BF_RB = bf_rb, p_upper = p_val[c('root')], mean = output$product[3], est_geno = output$product[4],
var_ran = output$product[5], rv_mean_es = re_rv_es, rv_sd_es = re_rv_sd, mean_cov=re_cov))


}
