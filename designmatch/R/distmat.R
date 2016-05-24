distmat <-
function (t_ind, X_mat, calip_cov = NULL, calip_size = NULL, calip_penalty = NULL, near_exact_covs = NULL, near_exact_penalties = NULL, digits = 1) {
    dist_mat = .smahal(t_ind, X_mat)
    dist_mat = dist_mat/mean(dist_mat)
    if (!is.null(calip_cov)) {
        dist_mat = .addcalip(dist_mat, t_ind, calip_cov, calip_size, calip_penalty)
    }
    if (!is.null(near_exact_covs)) {
    	for (i in 1:ncol(near_exact_covs)) {
    		penalty_mat = abs(outer(near_exact_covs[t_ind==1, i], near_exact_covs[t_ind==0, i], "-"))
    		penalty_mat = (penalty_mat!=0)*(near_exact_penalties[i])
    		dist_mat = dist_mat+penalty_mat
        }
    }
    dist_mat = round(dist_mat/mean(dist_mat), digits)
    dist_mat
}
