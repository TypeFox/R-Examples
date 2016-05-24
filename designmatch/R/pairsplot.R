pairsplot <-
function(cov1, cov2, t_id, c_id, xlab = "", ylab = "", main = "") {
    plot(cov1[t_id], cov2[t_id], xlim = c(min(cov1[c(t_id, c_id)]), max(cov1[c(t_id, c_id)])), ylim = c(min(cov2[c(t_id, c_id)]), max(cov2[c(t_id, c_id)])), xlab = xlab, ylab = ylab, main = main, pch = 1, col = "red")
    points(cov1[c_id], cov2[c_id], pch = 0, col = "blue")
    for (i in 1:length(t_id)) {
        segments(cov1[t_id[i]], cov2[t_id[i]], cov1[c_id[i]], 
				 cov2[c_id[i]], col = "grey")
    }
    abline(v = mean(cov1[t_id]), b = 0, col = "red")
    abline(v = mean(cov1[c_id]), b = 0, col = "blue")
    abline(h = mean(cov2[t_id]), b = 0, col = "red")
    abline(h = mean(cov2[c_id]), b = 0, col = "blue")
    legend_aux1 = c("Treated", "Controls")
    legend_aux2 = c("red", "blue")
    legend_aux3 = c(1, 0)
    legend("topright", legend_aux1, bty = "n", pch = legend_aux3, col = legend_aux2)
}
