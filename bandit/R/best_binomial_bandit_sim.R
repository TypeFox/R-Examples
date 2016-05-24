best_binomial_bandit_sim <- function(x,
 n,
 alpha = 1,
 beta = 1,
 ndraws = 5000
) {
    return(prob_winner(sim_post(x,n,alpha,beta,ndraws)))
}
