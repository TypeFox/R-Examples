mcmc_epsp <-
function(pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c, iter=10000, burn=500, sd=0.1)
{

num_m <- length(pheno_m)
num_d <- length(pheno_d)
num_a <- ncol(B_des_a_m)
num_c <- ncol(B_des_c_m)

B_a_m <- t(B_des_a_m)
B_a_d <- t(B_des_a_d)
B_c_m <- t(B_des_c_m)
B_c_d <- t(B_des_c_d)

if(var_b_a <= 0.05)
{
ei_a <- eigen(D_a)
G_a <- t(ei_a$vectors)
}else{
G_a <- diag(rep(1,num_a))
}

if(var_b_c <= 0.05)
{
ei_c <- eigen(D_c)
G_c <- t(ei_c$vectors)
}else{
G_c <- diag(rep(1,num_c))
}

multResult <- rep(0,num_a+num_c+(num_a+1)*num_a/2+(num_c+1)*num_c/2+1)

output =.C("CWrapper_mcmc",
product = as.double(multResult),
num_p_mz = as.integer(num_m),
num_p_dz = as.integer(num_d),
num_col_a = as.integer(num_a),
num_col_c = as.integer(num_c),
ph_m = as.double(pheno_m),
ph_d = as.double(pheno_d),
B_des_a_m = as.double(B_a_m),
B_des_a_d = as.double(B_a_d),
B_des_c_m = as.double(B_c_m),
B_des_c_d = as.double(B_c_d),
G_a = as.double(G_a),
G_c = as.double(G_c),
var = as.double(var),
var_b_a = as.double(var_b_a),
var_b_c = as.double(var_b_c),
D_a = as.integer(D_a),
D_c = as.integer(D_c),
iter_n = as.integer(iter),
burn = as.integer(burn),
sd_mcmc = as.double(sd)
)

beta_a_mc <- output$product[1:num_a]
beta_c_mc <- output$product[(1+num_a):(num_a+num_c)]

k <- 1
cov_a <- matrix(0, num_a, num_a)
for(i in 1:num_a)
{
for(j in i:num_a)
{
cov_a[i,j] <- output$product[num_a+num_c+k]
cov_a[j,i] <- cov_a[i,j]
k <- k + 1
}
}
cov_c <- matrix(0, num_c, num_c)
for(i in 1:num_c)
{
for(j in i:num_c)
{
cov_c[i,j] <- output$product[num_a+num_c+k]
cov_c[j,i] <- cov_c[i,j]
k <- k + 1
}
}


return(list(beta_a_mc = beta_a_mc, beta_c_mc = beta_c_mc, cov_a = cov_a, cov_c = cov_c))


}

mcmc_epsp_AtEt <-
function(pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, var_b_a, var_b_e, D_a, D_e, iter=10000, burn=500, sd=0.1)
{

num_m <- length(pheno_m)
num_d <- length(pheno_d)
num_a <- ncol(B_des_a_m)
num_e <- ncol(B_des_e_m)

B_a_m <- t(B_des_a_m)
B_a_d <- t(B_des_a_d)
B_e_m <- t(B_des_e_m)
B_e_d <- t(B_des_e_d)

if(var_b_a <= 0.05)
{
ei_a <- eigen(D_a)
#G_a <- diag(sqrt(c(ei_a$value[1:(num_a-2)],0,0)))%*%ei_a$vectors
#G_a[num_a-1,] <- rep(1,num_a)
#G_a[num_a,] <- 1:num_a
G_a <- t(ei_a$vectors)
}else{
G_a <- diag(rep(1,num_a))
}

if(var_b_e <= 0.05)
{
ei_e <- eigen(D_e)
#G_e <- diag(sqrt(c(ei_e$value[1:(num_e-2)],0,0)))%*%ei_e$vectors
#G_e[num_e-1,] <- rep(1,num_e)
#G_e[num_e,] <- 1:num_e
G_e <- t(ei_e$vectors)
}else{
G_e <- diag(rep(1,num_e))
}

multResult <- rep(0,num_a+num_e+(num_a+1)*num_a/2+(num_e+1)*num_e/2+1)
var <- -1

output =.C("CWrapper_mcmc",
product = as.double(multResult),
num_p_mz = as.integer(num_m),
num_p_dz = as.integer(num_d),
num_col_a = as.integer(num_a),
num_col_c = as.integer(num_e),
ph_m = as.double(pheno_m),
ph_d = as.double(pheno_d),
B_des_a_m = as.double(B_a_m),
B_des_a_d = as.double(B_a_d),
B_des_c_m = as.double(B_e_m),
B_des_c_d = as.double(B_e_d),
G_a = as.double(G_a),
G_c = as.double(G_e),
var = as.double(var),
var_b_a = as.double(var_b_a),
var_b_c = as.double(var_b_e),
D_a = as.integer(D_a),
D_c = as.integer(D_e),
iter_n = as.integer(iter),
burn = as.integer(burn),
sd_mcmc = as.double(sd)
)

beta_a_mc <- output$product[1:num_a]
beta_e_mc <- output$product[(1+num_a):(num_a+num_e)]

k <- 1
cov_a <- matrix(0, num_a, num_a)
for(i in 1:num_a)
{
for(j in i:num_a)
{
cov_a[i,j] <- output$product[num_a+num_e+k]
cov_a[j,i] <- cov_a[i,j]
k <- k + 1
}
}
cov_e <- matrix(0, num_e, num_e)
for(i in 1:num_e)
{
for(j in i:num_e)
{
cov_e[i,j] <- output$product[num_a+num_e+k]
cov_e[j,i] <- cov_e[i,j]
k <- k + 1
}
}

return(list(beta_a_mc = beta_a_mc, beta_e_mc = beta_e_mc, cov_a = cov_a, cov_e = cov_e))

}

mcmc_epsp_AtCtEt <-
function(pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, var_b_a, var_b_c, var_b_e, D_a, D_c, D_e, iter=10000, burn=500, sd=0.1)
{

num_m <- length(pheno_m)
num_d <- length(pheno_d)
num_a <- ncol(B_des_a_m)
num_c <- ncol(B_des_c_m)
num_e <- ncol(B_des_e_m)

B_a_m <- t(B_des_a_m)
B_a_d <- t(B_des_a_d)
B_c_m <- t(B_des_c_m)
B_c_d <- t(B_des_c_d)
B_e_m <- t(B_des_e_m)
B_e_d <- t(B_des_e_d)

if((var_b_a <= 0.05)&(num_a>2))
{
ei_a <- eigen(D_a)
G_a <- t(ei_a$vectors)
}else{
G_a <- diag(rep(1,num_a))
}

if((var_b_c <= 0.05)&(num_c>2))
{
ei_c <- eigen(D_c)
G_c <- t(ei_c$vectors)
}else{
G_c <- diag(rep(1,num_c))
}

if((var_b_e <= 0.05)&(num_e>2))
{
ei_e <- eigen(D_e)
G_e <- t(ei_e$vectors)
}else{
G_e <- diag(rep(1,num_e))
}

num_t <- num_a+num_c+num_e
# multResult <- rep(0,num_a+num_c+num_e+(num_a+1)*num_a/2+(num_c+1)*num_c/2+(num_e+1)*num_e/2+1)
multResult <- rep(0,num_t+(num_t+1)*num_t/2+1)

output =.C("CWrapper_mcmc_atctet",
product = as.double(multResult),
num_p_mz = as.integer(num_m),
num_p_dz = as.integer(num_d),
num_col_a = as.integer(num_a),
num_col_c = as.integer(num_c),
num_col_e = as.integer(num_e),
ph_m = as.double(pheno_m),
ph_d = as.double(pheno_d),
B_des_a_m = as.double(B_a_m),
B_des_a_d = as.double(B_a_d),
B_des_c_m = as.double(B_c_m),
B_des_c_d = as.double(B_c_d),
B_des_e_m = as.double(B_e_m),
B_des_e_d = as.double(B_e_d),
G_a = as.double(G_a),
G_c = as.double(G_c),
G_e = as.double(G_e),
var_b_a = as.double(var_b_a),
var_b_c = as.double(var_b_c),
var_b_e = as.double(var_b_e),
D_a = as.integer(D_a),
D_c = as.integer(D_c),
D_e = as.integer(D_e),
iter_n = as.integer(iter),
burn = as.integer(burn),
sd_mcmc = as.double(sd)
)

beta_a_mc <- output$product[1:num_a]
beta_c_mc <- output$product[(1+num_a):(num_a+num_c)]
beta_e_mc <- output$product[(1+num_a+num_c):(num_a+num_c+num_e)]

k <- 1
cov_t <- matrix(0, num_t, num_t)
for(i in 1:num_t)
{
for(j in i:num_t)
{
cov_t[i,j] <- output$product[num_t+k]
cov_t[j,i] <- cov_t[i,j]
k <- k + 1
}
}

return(list(beta_a_mc = beta_a_mc, beta_c_mc = beta_c_mc, beta_e_mc = beta_e_mc, cov = cov_t))


}
