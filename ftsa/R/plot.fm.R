`plot.fm` <- function (x, order, xlab1 = x$y$xname, ylab1 = "Principal component", 
    xlab2 = "Time", ylab2 = "Coefficient", mean.lab = "Mean", 
    level.lab = "Level", main.title = "Main effects", interaction.title = "Interaction", 
    basiscol = 1, coeffcol = 1, outlier.col = 2, outlier.pch = 19, 
    outlier.cex = 0.5, ...) 
{
       if(class(x)[1] == "ftsm"){
       oldpar <- par(no.readonly = TRUE)
       mean <- is.element("mean", colnames(x$basis))
       level <- is.element("level", colnames(x$basis))
       m <- mean + level
       order2 <- ncol(x$basis) - m
       if (missing(order)) 
           order <- order2
       n <- order + (mean | level)
       par(mfcol = c(2, n))
       if (mean) 
           plot(x$y$x, x$basis[, "mean"], type = "l", lty = 1, xlab = xlab1, 
                ylab = mean.lab, main = main.title, col = basiscol, 
                ...)
       if (level) 
           plot.ts(x$coeff[, "level"], xlab = xlab2, ylab = level.lab, 
                 main = ifelse(mean, "", main.title), col = coeffcol, 
                 ...)
       if (m == 1) 
           plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
                xlab = "", ylab = "")
       if (order > 0) {
           for (i in 1:order) {
                yl1 <- ifelse(n > 1, paste(ylab1, i), ylab1)
                yl2 <- ifelse(n > 1, paste(ylab2, i), ylab2)
                plot(x$y$x, x$basis[, m + i], type = "l", lty = 1, 
                     xlab = xlab1, ylab = yl1, col = basiscol, ...)
                if (i == 1) 
                    title(interaction.title)
                    plot.ts(x$coeff[, m + i], xlab = xlab2, ylab = yl2, 
                        col = coeffcol, ...)
                if (!is.null(x$weights)) 
                     points(time(x$coeff)[x$weights < 0.1], x$coeff[x$weights < 
                   0.1, i + m], pch = outlier.pch, col = outlier.col, 
                  cex = outlier.cex)
           }
        }
        par(oldpar)
      }
      if(class(x)[1] == "fm"){
         oldpar <- par(no.readonly = TRUE)
         if(missing(order)){
            order = dim(x$T)[2]
         }
         sco = x$T
         basis = x$Q
         par(mfrow = c(2, order + 1))
         plot(x$meanY, ylab = mean.lab, main = main.title)
         plot(basis[,1], xlab = xlab1, ylab = paste(ylab1, 1), type = "l", main = interaction.title)          
         for(i in 1:(order-1)){
             plot(basis[,i+1], xlab = xlab1, ylab = paste(ylab1, i+1), type = "l")
         }  
         plot(basis[,1], type = "n", xlab = "", ylab = "", axes = FALSE)     
         for(i in 1:order){
             plot(ts(sco[,i], start = tsp(x$y$time)[1], end = tsp(x$y$time)[2]), xlab = xlab2, ylab = paste(ylab2, i), type = "l")
         }
         par(oldpar)
      } 
}
