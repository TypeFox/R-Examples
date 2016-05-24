loveplot <-
function (X_mat, t_id, c_id, v_line, legend_position = "topright") {
	X_mat_t = X_mat[t_id, ]
	X_mat_c_before = X_mat[-t_id, ]
	X_mat_c_before_mean = apply(X_mat_c_before, 2, mean)
	X_mat_t_mean = apply(X_mat_t, 2, mean)
	X_mat_t_var = apply(X_mat_t, 2, var)
	X_mat_c_before_var = apply(X_mat_c_before, 2, var)
	std_dif_before = (X_mat_t_mean - X_mat_c_before_mean)/sqrt((X_mat_t_var + X_mat_c_before_var)/2)
	X_mat_c_after = X_mat[c_id, ]
	X_mat_c_after_mean = apply(X_mat_c_after, 2, mean)
	std_dif_after = (X_mat_t_mean - X_mat_c_after_mean)/sqrt((X_mat_t_var + X_mat_c_before_var)/2)
	#library("lattice")
	abs_std_dif_before = abs(std_dif_before)
	n_aux = length(abs_std_dif_before)
	abs_std_dif_after = abs(std_dif_after)
	dotchart(abs_std_dif_before[n_aux:1], labels = colnames(X_mat)[n_aux:1], cex = 0.7, pch = "", color = , main = "", xlim = c(0, 1), xlab = "Absolute standardized differences in means")
	points(abs_std_dif_before[n_aux:1], y = 1:ncol(X_mat), cex = 0.9, pch = 0)
	points(abs_std_dif_after[n_aux:1], y = 1:ncol(X_mat), cex = 0.8, pch = 8, col = "blue")
	legend(legend_position, c("Before matching", "After matching"), cex = 0.8, bty = "n", pch = c(0, 8), col = c("black", "blue"))
	abline(v = v_line, lty = 2)
}
