SVDplot = function (object, order = 3, plot = TRUE, plot.type = c("fts", 
    "image"), mfrow = c(2, 3)) 
{
    data = object$y
    p = dim(data)[1]
    n = dim(data)[2]
    pcscore = svd(t(data))$u[, 1:order]
    pc = svd(t(data))$v[, 1:order]
    d = svd(t(data))$d[1:order]
    svdr = array(, dim = c(p, n, order))
    if (plot == TRUE) {
        plot.type = match.arg(plot.type)
        xlab = object$xname
        par(mfrow = mfrow)
        if (plot.type == "fts") {
            plot(object, main = "Original data")
            for (i in 1:order) {
                svdr[, , i] = as.matrix(pc[, i] * d[i]) %*% t(pcscore[, i])
                plot(fts(object$x, svdr[, , i]), ylab = "", xlab = xlab, 
                  main = paste(expression(SVD), i, sep = ""))
            }
            svdrecon = apply(svdr, c(1, 2), sum)
            resi = data - svdrecon
            plot(fts(object$x, svdrecon), ylab = "", xlab = xlab, 
                main = "Reconstruction")
            plot(fts(object$x, resi), ylab = "", xlab = xlab, 
                main = "Residual")
        }
        else {
            image(object$x, as.numeric(colnames(data)), data, 
                main = "Original data", xlab = xlab, ylab = object$yname)
            box()
            for (i in 1:order) {
                svdr[, , i] = as.matrix(pc[, i] * d[i]) %*% t(pcscore[, 
                  i])
                image(object$x, as.numeric(colnames(data)), svdr[, 
                  , i], ylab = "", xlab = xlab, main = paste(expression(SVD), 
                  i, sep = ""))
                box()
            }
            svdrecon = apply(svdr, c(1, 2), sum)
            resi = data - svdrecon
            image(object$x, as.numeric(colnames(data)), svdrecon, 
                ylab = "", xlab = xlab, main = "Reconstruction")
            box()
            image(object$x, as.numeric(colnames(data)), resi, 
                ylab = "", xlab = xlab, main = "Residual")
            box()
        }
    }
    else {
		svdr = array(, dim=c(p,n,order))
         for (i in 1:order) 
	   {
              svdr[,,i] = as.matrix(pc[, i] * d[i]) %*% t(pcscore[, i])
	   }	
	   svdrecon = apply(svdr, c(1, 2), sum)
         resi = data - svdrecon
         return(list(svds = svdr, reconstruction = svdrecon, residual = resi))
    }
}
