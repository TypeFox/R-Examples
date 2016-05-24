ecdfplot <-
function(x, t_id, c_id, main_title = "", legend_position = "right") {
	xlim_min = quantile(x, .01)
	xlim_max = quantile(x, .99)
	aux1 = ecdf(x[t_id])
	aux2 = ecdf(x[c_id])
	if (length(c_id) > 250) {
		aux2 = ecdf(sample(x, 250))
	}
	plot(aux1, xlim = c(xlim_min, xlim_max), ylim = c(0, 1), xlab = "x", ylab = "CDF(x)", main = main_title, col = "black", pch = 0)
	par(new=TRUE)
	plot(aux2, xlim = c(xlim_min, xlim_max), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", main = "", col = "blue", pch = 8)
	legend_aux1 = c("Treated", "Controls")
	legend_aux2 = c("black", "blue")
	legend_aux3 = c(0, 8)
	legend(legend_position, legend_aux1, bty = "n", pch = legend_aux3, col = legend_aux2)
	par(new=FALSE)	
}
