
# Adding 40% Uncertainty to fixed effects log-normal (not Favail)
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                      bpop_vals,
                      ones(length(bpop_vals),1)*(bpop_vals*0.4)^2) # 40% of bpop value
bpop_vals_ed_ln["Favail",]  <- c(0,1,0)

pars.ln <- pargen(par=bpop_vals_ed_ln,
               user_dist_pointer=NULL,
               sample_size=1000,
               bLHS=1,
               sample_number=NULL,
               poped.db)


# Adding 10% Uncertainty to fixed effects normal-distribution (not Favail)
bpop_vals_ed_n <- cbind(ones(length(bpop_vals),1)*1, # log-normal distribution
                      bpop_vals,
                      ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed_n["Favail",]  <- c(0,1,0)

pars.n <- pargen(par=bpop_vals_ed_n,
               user_dist_pointer=NULL,
               sample_size=1000,
               bLHS=1,
               sample_number=NULL,
               poped.db)


# Adding 10% Uncertainty to fixed effects uniform-distribution (not Favail)
bpop_vals_ed_u <- cbind(ones(length(bpop_vals),1)*2, # uniform distribution
                        bpop_vals,
                        ones(length(bpop_vals),1)*(bpop_vals*0.1)) # 10% of bpop value
bpop_vals_ed_u["Favail",]  <- c(0,1,0)

pars.u <- pargen(par=bpop_vals_ed_u,
                 user_dist_pointer=NULL,
                 sample_size=1000,
                 bLHS=1,
                 sample_number=NULL,
                 poped.db)


# Adding user defined distributions
bpop_vals_ed_ud <- cbind(ones(length(bpop_vals),1)*3, # user dfined distribution
                         bpop_vals,
                         bpop_vals*0.1) # 10% of bpop value
bpop_vals_ed_ud["Favail",]  <- c(0,1,0)

# A normal distribution
my_dist <- function(...){
  par_vec <- rnorm(c(1,1,1,1),mean=bpop_vals_ed_ud[,2],sd=bpop_vals_ed_ud[,3])
}

pars.ud <- pargen(par=bpop_vals_ed_ud,
                  user_dist_pointer=my_dist,
                  sample_size=1000,
                  bLHS=1,
                  sample_number=NULL,
                  poped.db)


