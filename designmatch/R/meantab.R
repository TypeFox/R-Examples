meantab <-
function(X_mat, t_ind, t_id, c_id, exact = NULL, digits = 2) {
#t_ind = c(rep(1, length(t_id)), rep(0, length(c_id)))
#X_mat = X_mat[c(t_id, c_id), ]
	n_vrs = ncol(X_mat)
	n_obs = nrow(X_mat)
	output = matrix(NA, n_vrs, 7)
	rownames(output) = colnames(X_mat)
	colnames(output) = c("Mis", "Min", "Max", "Mean T", "Mean C", "Std Dif", "P-val")
	for (j in 1:n_vrs) {
		n_miss = sum(is.na(X_mat[, j]))
		miss_pct = round(n_miss/n_obs, 2)
		if (n_miss < n_obs) {
#yes = unlist(X_mat[t_ind==1, j])
#no = unlist(X_mat[t_ind==0, j])
#mean_yes = mean(yes, na.rm = TRUE)
#mean_no = mean(no, na.rm = TRUE)
#min_val = min(c(yes, no), na.rm = TRUE)
#max_val = max(c(yes, no), na.rm = TRUE)
#pooled_sd = sqrt((var(yes, na.rm = TRUE)+var(no, na.rm = TRUE))/2)
			yes = unlist(X_mat[t_id, j])
			no = unlist(X_mat[c_id, j])
			
			yes_before = unlist(X_mat[t_ind==1, j])
			no_before = unlist(X_mat[t_ind==0, j])
			
			mean_yes = mean(yes, na.rm = TRUE)
			mean_no = mean(no, na.rm = TRUE)
			min_val = min(c(yes, no), na.rm = TRUE)
			max_val = max(c(yes, no), na.rm = TRUE)
			pooled_sd = sqrt((var(yes_before, na.rm = TRUE)+var(no_before, na.rm = TRUE))/2)
			std_diff = NA
			p_val = NA
			if (max_val > min_val) {
				std_diff = (mean_yes-mean_no)/pooled_sd
				if (is.null(exact)) {
					p_val = t.test(yes, no)$p.value
				}
				if (!is.null(exact)) {
					if (exact[j]=="w") {
						p_val = wilcox.test(yes, no)$p.value
					}
					if (exact[j]=="f") {
						tab_aux = table(t_ind[c(t_id, c_id)], unlist(X_mat[c(t_id, c_id), j]) )
						if (ncol(tab_aux)==1) {
							p_val = 1
						}
						if (ncol(tab_aux)!=1) {
							p_val = fisher.test(tab_aux)$p.value
						}
					}	
				}
			}
			output[j, ] = c(miss_pct, min_val, max_val, mean_yes, mean_no, std_diff, p_val)
		}
	}
	round(output, digits)
}
