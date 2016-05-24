loglik_AtCtEt_esp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d)
{

beta_a <- param[1:ncol(B_des_a_m)]
beta_c <- param[(1+ncol(B_des_a_m)):(ncol(B_des_c_m)+ncol(B_des_a_m))]
beta_e <- param[(1+ncol(B_des_a_m)+ncol(B_des_c_m)):(ncol(B_des_e_m)+ncol(B_des_c_m)+ncol(B_des_a_m))]

nll <- .Call('loglik_AtCtEt_esp_c', beta_a, beta_c, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d)

return(nll)
}

gr_AtCtEt_esp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d)
{
beta_a <- param[1:ncol(B_des_a_m)]
beta_c <- param[(1+ncol(B_des_a_m)):(ncol(B_des_c_m)+ncol(B_des_a_m))]
beta_e <- param[(1+ncol(B_des_a_m)+ncol(B_des_c_m)):(ncol(B_des_e_m)+ncol(B_des_c_m)+ncol(B_des_a_m))]

num_beta_c <- ncol(B_des_c_m)
num_beta_a <- ncol(B_des_a_m)
num_beta_e <- ncol(B_des_e_m)

d <- .Call('gr_AtCtEt_esp_c', beta_a, beta_c, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d)

return(d) 
}
