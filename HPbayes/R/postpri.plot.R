postpri.plot <-
function(..., prior, hpp, box = FALSE, type = "l", line.col = c("black", 
    "red"), line.bound = TRUE, rowcol = c(2, 4)) 
{
    par(mfrow = rowcol, ...)
    lc <- as.matrix(line.col)
    type1 <- type
    type2 <- type
    pri.a <- density(prior[, 1])
    pri.ax <- pri.a$x
    pri.ay <- pri.a$y
    post.a <- density(hpp[, 1])
    post.ax <- post.a$x
    post.ay <- post.a$y
    pri.b <- density(prior[, 2])
    pri.bx <- pri.b$x
    pri.by <- pri.b$y
    post.b <- density(hpp[, 2])
    post.bx <- post.b$x
    post.by <- post.b$y
    pri.c <- density(prior[, 3])
    pri.cx <- pri.c$x
    pri.cy <- pri.c$y
    post.c <- density(hpp[, 3])
    post.cx <- post.c$x
    post.cy <- post.c$y
    pri.d <- density(prior[, 4])
    pri.dx <- pri.d$x
    pri.dy <- pri.d$y
    post.d <- density(hpp[, 4])
    post.dx <- post.d$x
    post.dy <- post.d$y
    pri.e <- density(prior[, 5])
    pri.ex <- pri.e$x
    pri.ey <- pri.e$y
    post.e <- density(hpp[, 5])
    post.ex <- post.e$x
    post.ey <- post.e$y
    pri.f <- density(prior[, 6])
    pri.fx <- pri.f$x
    pri.fy <- pri.f$y
    post.f <- density(hpp[, 6])
    post.fx <- post.f$x
    post.fy <- post.f$y
    pri.g <- density(prior[, 7])
    pri.gx <- pri.g$x
    pri.gy <- pri.g$y
    post.g <- density(hpp[, 7])
    post.gx <- post.g$x
    post.gy <- post.g$y
    pri.h <- density(prior[, 8])
    pri.hx <- pri.h$x
    pri.hy <- pri.h$y
    post.h <- density(hpp[, 8])
    post.hx <- post.h$x
    post.hy <- post.h$y
    if (box == FALSE & line.bound == FALSE) {
        plot(pri.ax, pri.ay, type = type1, ylim = c(min(pri.ay, 
            post.ay), max(pri.ay, post.ay)), xlim = c(min(pri.ax, 
            post.ax), max(pri.ax, post.ax)), ylab = "density", 
            xlab = "Prior/Posterior", main = "A", col = lc[1, 
                ])
        points(post.ax, post.ay, type = type2, col = lc[2, ])
        plot(pri.bx, pri.by, type = type1, ylim = c(min(pri.by, 
            post.by), max(pri.by, post.by)), xlim = c(min(pri.bx, 
            post.bx), max(pri.bx, post.bx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "B", col = lc[1, 
                ])
        points(post.bx, post.by, type = type2, col = lc[2, ])
        plot(pri.cx, pri.cy, type = type1, ylim = c(min(pri.cy, 
            post.cy), max(pri.cy, post.cy)), xlim = c(min(pri.cx, 
            post.cx), max(pri.cx, post.cx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "C", col = lc[1, 
                ])
        points(post.cx, post.cy, type = type2, col = lc[2, ])
        plot(pri.dx, pri.dy, type = type1, ylim = c(min(pri.dy, 
            post.dy), max(pri.dy, post.dy)), xlim = c(min(pri.dx, 
            post.dx), max(pri.dx, post.dx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "D", col = lc[1, 
                ])
        points(post.dx, post.dy, type = type2, col = lc[2, ])
        plot(pri.ex, pri.ey, type = type1, ylim = c(min(pri.ey, 
            post.ey), max(pri.ey, post.ey)), xlim = c(min(pri.ex, 
            post.ex), max(pri.ex, post.ex)), ylab = "density", 
            xlab = "Prior/Posterior", main = "E", col = lc[1, 
                ])
        points(post.ex, post.ey, type = type2, col = lc[2, ])
        plot(pri.fx, pri.fy, type = type1, ylim = c(min(pri.fy, 
            post.fy), max(pri.fy, post.fy)), xlim = c(min(pri.fx, 
            post.fx), max(pri.fx, post.fx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "F", col = lc[1, 
                ])
        points(post.fx, post.fy, type = type2, col = lc[2, ])
        plot(pri.gx, pri.gy, type = type1, ylim = c(min(pri.gy, 
            post.gy), max(pri.gy, post.gy)), xlim = c(min(pri.gx, 
            post.gx), max(pri.gx, post.gx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "G", col = lc[1, 
                ])
        points(post.gx, post.gy, type = type2, col = lc[2, ])
        plot(pri.hx, pri.hy, type = type1, ylim = c(min(pri.hy, 
            post.hy), max(pri.hy, post.hy)), xlim = c(min(pri.hx, 
            post.hx), max(pri.hx, post.hx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "H", col = lc[1, 
                ])
        points(post.hx, post.hy, type = type2, col = lc[2, ])
    }
    
    if(box == FALSE & line.bound == TRUE){
    	plot(post.ax, post.ay, type = type2, ylim = c(min(pri.ay, 
            post.ay), max(pri.ay, post.ay)), xlim = c(min(pri.ax, 
            post.ax), max(pri.ax, post.ax)), ylab = "density", 
            xlab = "Prior/Posterior", main = "A", col = lc[1,])
        segments(x0=min(prior[,1]), y0=min(pri.ay), x1=min(prior[,1]), y1=max(pri.ay), lty=2, col = lc[2,])
        segments(x0=min(prior[,1]), y0=max(pri.ay), x1=max(prior[,1]), y1=max(pri.ay), lty=2, col = lc[2,])
        segments(x0=max(prior[,1]), y0=max(pri.ay), x1=max(prior[,1]), y1=min(pri.ay), lty=2, col = lc[2,])
        
        plot(post.bx, post.by, type = type2, ylim = c(min(pri.by, 
            post.by), max(pri.by, post.by)), xlim = c(min(pri.bx, 
            post.bx), max(pri.bx, post.bx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "B", col = lc[1,])
        segments(x0=min(prior[,2]), y0=min(pri.by), x1=min(prior[,2]), y1=max(pri.by), lty=2, col = lc[2,])
        segments(x0=min(prior[,2]), y0=max(pri.by), x1=max(prior[,2]), y1=max(pri.by), lty=2, col = lc[2,])
        segments(x0=max(prior[,2]), y0=max(pri.by), x1=max(prior[,2]), y1=min(pri.by), lty=2, col = lc[2,])
        
		plot(post.cx, post.cy, type = type2, ylim = c(min(pri.cy, 
            post.cy), max(pri.cy, post.cy)), xlim = c(min(pri.cx, 
            post.cx), max(pri.cx, post.cx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "C", col = lc[1,])
        segments(x0=min(prior[,3]), y0=min(pri.cy), x1=min(prior[,3]), y1=max(pri.cy), lty=2, col = lc[2,])
        segments(x0=min(prior[,3]), y0=max(pri.cy), x1=max(prior[,3]), y1=max(pri.cy), lty=2, col = lc[2,])
        segments(x0=max(prior[,3]), y0=max(pri.cy), x1=max(prior[,3]), y1=min(pri.cy), lty=2, col = lc[2,])
                
        plot(post.dx, post.dy, type = type2, ylim = c(min(pri.dy, 
            post.dy), max(pri.dy, post.dy)), xlim = c(min(pri.dx, 
            post.dx), max(pri.dx, post.dx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "D", col = lc[1,])
        segments(x0=min(prior[,4]), y0=min(pri.dy), x1=min(prior[,4]), y1=max(pri.dy), lty=2, col = lc[2,])
        segments(x0=min(prior[,4]), y0=max(pri.dy), x1=max(prior[,4]), y1=max(pri.dy), lty=2, col = lc[2,])
        segments(x0=max(prior[,4]), y0=max(pri.dy), x1=max(prior[,4]), y1=min(pri.dy), lty=2, col = lc[2,])
                
        plot(post.ex, post.ey, type = type2, ylim = c(min(pri.ey, 
            post.ey), max(pri.ey, post.ey)), xlim = c(min(pri.ex, 
            post.ex), max(pri.ex, post.ex)), ylab = "density", 
            xlab = "Prior/Posterior", main = "E", col = lc[1,])
        segments(x0=min(prior[,5]), y0=min(pri.ey), x1=min(prior[,5]), y1=max(pri.ey), lty=2, col = lc[2,])
        segments(x0=min(prior[,5]), y0=max(pri.ey), x1=max(prior[,5]), y1=max(pri.ey), lty=2, col = lc[2,])
        segments(x0=max(prior[,5]), y0=max(pri.ey), x1=max(prior[,5]), y1=min(pri.ey), lty=2, col = lc[2,])
                
        plot(post.fx, post.fy, type = type2, ylim = c(min(pri.fy, 
            post.fy), max(pri.fy, post.fy)), xlim = c(min(pri.fx, 
            post.fx), max(pri.fx, post.fx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "F", col = lc[1,])
        segments(x0=min(prior[,6]), y0=min(pri.fy), x1=min(prior[,6]), y1=max(pri.fy), lty=2, col = lc[2,])
        segments(x0=min(prior[,6]), y0=max(pri.fy), x1=max(prior[,6]), y1=max(pri.fy), lty=2, col = lc[2,])
        segments(x0=max(prior[,6]), y0=max(pri.fy), x1=max(prior[,6]), y1=min(pri.fy), lty=2, col = lc[2,])
                
        plot(post.gx, post.gy, type = type2, ylim = c(min(pri.gy, 
            post.gy), max(pri.gy, post.gy)), xlim = c(min(pri.gx, 
            post.gx), max(pri.gx, post.gx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "G", col = lc[1,])
        segments(x0=min(prior[,7]), y0=min(pri.gy), x1=min(prior[,7]), y1=max(pri.gy), lty=2, col = lc[2,])
        segments(x0=min(prior[,7]), y0=max(pri.gy), x1=max(prior[,7]), y1=max(pri.gy), lty=2, col = lc[2,])
        segments(x0=max(prior[,7]), y0=max(pri.gy), x1=max(prior[,7]), y1=min(pri.gy), lty=2, col = lc[2,])
                
        plot(post.hx, post.hy, type = type2, ylim = c(min(pri.hy, 
            post.hy), max(pri.hy, post.hy)), xlim = c(min(pri.hx, 
            post.hx), max(pri.hx, post.hx)), ylab = "density", 
            xlab = "Prior/Posterior", main = "H", col = lc[1,])
        segments(x0=min(prior[,8]), y0=min(pri.hy), x1=min(prior[,8]), y1=max(pri.hy), lty=2, col = lc[2,])
        segments(x0=min(prior[,8]), y0=max(pri.hy), x1=max(prior[,8]), y1=max(pri.hy), lty=2, col = lc[2,])
        segments(x0=max(prior[,8]), y0=max(pri.hy), x1=max(prior[,8]), y1=min(pri.hy), lty=2, col = lc[2,])   	}
    
    if (box == TRUE) {
        boxplot(prior[, 1], hpp[, 1], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "A", boxwex = 0.8)
        boxplot(prior[, 2], hpp[, 2], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "B", boxwex = 0.8)
        boxplot(prior[, 3], hpp[, 3], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "C", boxwex = 0.8)
        boxplot(prior[, 4], hpp[, 4], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "D", boxwex = 0.8)
        boxplot(prior[, 5], hpp[, 5], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "E", boxwex = 0.8)
        boxplot(prior[, 6], hpp[, 6], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "F", boxwex = 0.8)
        boxplot(prior[, 7], hpp[, 7], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "G", boxwex = 0.8)
        boxplot(prior[, 8], hpp[, 8], col = lc[1, ], names = c("Prior", 
            "Posterior"), main = "H", boxwex = 0.8)
    }
}

