plot.npmr <-
function(x, lambda, feature.names = TRUE, ...) {

    l = which.min((x$lambda - lambda)^2)
    B = x$B[,,l]
    SVD = svd(B)

    if (is.null(rownames(B))) rownames(B) = 1:nrow(B)
    if (is.null(colnames(B))) colnames(B) = 1:ncol(B)

    plot(SVD$u[, 1], SVD$u[, 2], col = 'darkorange', axes = FALSE, type = 'n',
        xlab = 'First latent feature', ylab = 'Second latent feature', 
        xlim = range(c(SVD$u[, 1], SVD$v[, 1])),
        ylim = range(c(SVD$u[, 2], SVD$v[, 2])), ...)
    if (feature.names) {
        text(SVD$u[, 1], SVD$u[, 2], labels = rownames(B), col = 'darkorange')
    } else {
        points(SVD$u[, 1], SVD$u[, 2], col = 'darkorange')
    }
    arrows(0, 0, SVD$v[, 1], SVD$v[, 2], length = 0.1, lty = 2, lwd = 0.5,
        col = 'dodgerblue')
    text(SVD$v[, 1], SVD$v[, 2], labels = colnames(B), col = 'dodgerblue',
        cex = 1.5)
    legend('bottomleft', c('Features', 'Classes'),
        text.col = c('darkorange', 'dodgerblue'))
}
