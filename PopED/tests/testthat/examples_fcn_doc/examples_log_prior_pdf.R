
# Adding 40% Uncertainty to fixed effects log-normal (not Favail)
bpop_vals <- c(CL=0.15, V=8, KA=1.0, Favail=1)
bpop_vals_ed_ln <- cbind(ones(length(bpop_vals),1)*4, # log-normal distribution
                         bpop_vals,
                         ones(length(bpop_vals),1)*(bpop_vals*0.4)^2) # 40% of bpop value
bpop_vals_ed_ln["Favail",]  <- c(0,1,0)
bpop_vals_ed_ln

# then compute the log density 
alpha <- bpop_vals_ed_ln[bpop_vals_ed_ln[,1]!=0,2]
log_prior_pdf(alpha=alpha, bpopdescr=bpop_vals_ed_ln, ddescr=poped.db$parameters$d)


# Adding 10% Uncertainty to fixed effects normal-distribution (not Favail)
bpop_vals_ed_n <- cbind(ones(length(bpop_vals),1)*1, # log-normal distribution
                        bpop_vals,
                        ones(length(bpop_vals),1)*(bpop_vals*0.1)^2) # 10% of bpop value
bpop_vals_ed_n["Favail",]  <- c(0,1,0)
bpop_vals_ed_n

# then compute the log density from log_prior_pdf
alpha <- bpop_vals_ed_n[bpop_vals_ed_n[,1]!=0,2]
log_prior_pdf(alpha=alpha, bpopdescr=bpop_vals_ed_n, ddescr=poped.db$parameters$d)


# Adding 10% Uncertainty to fixed effects uniform-distribution (not Favail)
bpop_vals_ed_u <- cbind(ones(length(bpop_vals),1)*2, # uniform distribution
                        bpop_vals,
                        ones(length(bpop_vals),1)*(bpop_vals*0.1)) # 10% of bpop value
bpop_vals_ed_u["Favail",]  <- c(0,1,0)
bpop_vals_ed_u

# then compute the log density 
alpha <- bpop_vals_ed_ln[bpop_vals_ed_u[,1]!=0,2]
log_prior_pdf(alpha=alpha, bpopdescr=bpop_vals_ed_u, ddescr=poped.db$parameters$d)

