loglik_AtCtEt_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c, B_des_e_m, B_des_e_d, beta_e, D_e)
{
var_b_a <- param[1]
var_b_c <- param[2]
var_b_e <- param[3]

nll <- .Call('loglik_AtCtEt_epsp_c', var_b_a, var_b_c, var_b_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c, B_des_e_m, B_des_e_d, beta_e, D_e)

return(nll)
}

loglik_AtCtEt_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, var_b_a, var_b_c, var_b_e, D_a, D_c, D_e)
{
beta_a <- param[1:ncol(D_a)]
beta_c <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_c))]
beta_e <- param[(ncol(D_a)+ncol(D_c)+1):(ncol(D_a)+ncol(D_c)+ncol(D_e))]

nll <- .Call('loglik_AtCtEt_epsp_g_c', beta_a, beta_c, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, var_b_a, var_b_c, var_b_e, D_a, D_c, D_e)

return(nll)
}

gr_AtCtEt_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c, B_des_e_m, B_des_e_d, beta_e, D_e)
{
var_b_a <- param[1]
var_b_c <- param[2]
var_b_e <- param[3]

d <- .Call('gr_AtCtEt_epsp_c', var_b_a, var_b_c, var_b_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c, B_des_e_m, B_des_e_d, beta_e, D_e)

return(d)
}

gr_AtCtEt_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, var_b_a, var_b_c, var_b_e, D_a, D_c, D_e)
{
beta_a <- param[1:ncol(D_a)]
beta_c <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_c))]
beta_e <- param[(ncol(D_a)+ncol(D_c)+1):(ncol(D_a)+ncol(D_c)+ncol(D_e))]

d <- .Call('gr_AtCtEt_epsp_g_c', beta_a, beta_c, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, var_b_a, var_b_c, var_b_e, D_a, D_c, D_e)

return(d)
}
