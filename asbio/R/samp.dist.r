dirty.dist <- function(s.size, parent = expression(rnorm(1)), cont = expression(rnorm(1, mean = 10)), prop.cont = 0.1){
vec <- seq(1,s.size)
for(i in 1 : s.size){
pure <- eval(parent)
dirty <- eval(cont)
vec[i]<- sample(c(pure, dirty), size = 1, prob = c(1 - prop.cont, prop.cont))
}
vec
} 

samp.dist<- function (parent = NULL, parent2 = NULL, biv.parent = NULL, s.size = 1, s.size2 = NULL, R = 1000, nbreaks = 50, stat = mean, stat2 = NULL, 
stat3 = NULL, stat4 = NULL, xlab = expression(bar(x)), func = NULL, show.n = TRUE, show.SE = FALSE, anim = TRUE, interval = 0.01, col.anim = "rainbow", digits = 3, ...) 
{
 
    if (is.null(col.anim)) 
        clr = NULL
    if (!is.null(col.anim)) {
        if (col.anim == "rainbow") 
            clr <- rainbow(R)
        if (col.anim == "heat.colors") {
            c1 <- heat.colors(R, 1)
            clr <- c1[length(c1):1]
        }
        if (col.anim == "gray") 
            clr <- gray(1 - (1:R)/R)
        if (col.anim != "rainbow" & col.anim != "gray" & col.anim != 
            "heat.colors") 
            clr <- rep(col.anim, R)
    }

if(!is.null(biv.parent)){
     
        func.res <- seq(1:R)
        for (i in 1:R) {
            if(is.expression(biv.parent)) part <- eval(biv.parent)
            if(!is.expression(biv.parent)){ 
                part <- sample(1:nrow(biv.parent), s.size, replace = FALSE)
                samp <- sample(1:nrow(biv.parent), s.size, replace = FALSE)
                part <- biv.parent[samp,]}
            
            func.res[i] <- func(part[,1],part[,2])
}
}
    
if(is.null(biv.parent)){    
    if (is.null(stat2) & is.null(stat3)) {
        func.res <- matrix(ncol = 1, nrow = R)
            for (i in 1:R) {
                if(is.expression(parent))func.res[i] <- stat(eval(parent))
                if(!is.expression(parent)){
                    part <- sample(parent, size = s.size, replace = FALSE)
                    func.res[i] <- stat(part)}
           }
        if (!is.null(func)) 
            func.res <- func(func.res, s.dist2 = NULL, s.size, s.size2 = NULL)
    }

    if (!is.null(stat2) | !is.null(stat3)){
        s.dist1 <- matrix(ncol = 1, nrow = R)
        s.dist2 <- matrix(ncol = 1, nrow = R)
        s.dist3 <- matrix(ncol = 1, nrow = R)
        s.dist4 <- matrix(ncol = 1, nrow = R)
         for(i in 1 : R){
             if(is.expression(parent)){
                sample1 <- eval(parent)
                if(!is.null(parent2)) sample2 <- eval(parent2)
                }
             if(!is.expression(parent)){
                sample1 <- sample(parent, size = s.size, replace = FALSE)
                if(!is.null(parent2)) sample2 <- sample(parent2, size = s.size, replace = FALSE)
                }
                        
                    s.dist1[i] <- stat(sample1)
                if (!is.null(stat2)) 
                    s.dist2[i] <- stat2(sample2)
                if (!is.null(stat3)) 
                    s.dist3[i] <- stat3(sample1)
                if (!is.null(stat4)) 
                    s.dist4[i] <- stat4(sample2)
              } 
        
        if (!is.null(stat2) & (is.null(stat3) | is.null(stat4))) {
            func.res <- func(s.dist1, s.dist2, s.size, s.size2)
        }
        if (!is.null(stat3) & (is.null(stat2) | is.null(stat4))) {
            func.res <- func(s.dist1, s.dist3, s.size, s.size2)
        }                                                               
        if (!is.null(stat2) & !is.null(stat3) & !is.null(stat4)) {
            func.res <- func(s.dist1, s.dist2, s.dist3, s.dist4, s.size, s.size2)
        }
        }
        }
        SE <- round(ifelse(is.matrix(func.res),apply(func.res, 2, sd),sd(func.res)),digits)
        brks <- seq(min(func.res),max(func.res),length.out=nbreaks)
        
        if (anim == TRUE) {
            for (i in 1:R) {
                dev.hold()
                hist(func.res, xlab = xlab,  
                  main = "", freq = FALSE, breaks = brks, border = "white",...)
                points(suppressWarnings(hist(func.res[1:i], plot = FALSE, 
                  breaks = brks, freq = FALSE)$mids), suppressWarnings(hist(func.res[1:i], 
                  plot = FALSE, breaks = brks, freq = FALSE)$density), 
                  type = "h", col = clr[i], lwd = 5)
                if(i < R) legend("topright",legend=i, bty="n", title="Sample")
                dev.flush()
                Sys.sleep(interval)
            }
        }
        if (anim == FALSE) {
            hist(func.res, xlab = xlab,  
                main = "", freq = FALSE, breaks = brks, ...)
        }
    
    if (show.SE == FALSE & show.n == TRUE) {
       if(is.null(s.size2))legend("topright", legend = paste("n = ", s.size), bty = "n")
       if(!is.null(s.size2))legend("topright", legend = c(paste("n1 = ", s.size), paste("n2 = ", s.size2)), bty = "n")
    }
    if (show.SE == TRUE & show.n == FALSE) {
       legend("topright", legend = paste("SE = ", SE), bty = "n")
    }
    if (show.SE == TRUE & show.n == TRUE) {
       if(is.null(s.size2)){
       legend("topright", legend = c(paste("n = ", s.size), 
            paste("SE = ", SE)), bty = "n")}
       if(!is.null(s.size2)){
       legend("topright", legend = c(paste("n1 = ", s.size), paste("n2 = ", s.size2), 
            paste("SE = ", SE)), bty = "n")}     
    }
}

#----------------------------- changing n ------------------------------------#

samp.dist.n<-function (parent, R = 500, n.seq = seq(1, 30), stat = mean, xlab = expression(bar(x)), 
    nbreaks = 50, func = NULL, show.n = TRUE, 
    show.SE = FALSE, est.density = TRUE, col.density = 4, lwd.density = 2, 
    est.ylim = TRUE, ylim = NULL, anim = TRUE, interval = 0.5, 
    col.anim = NULL, digits = 3, ...) 
{
    n.max <- max(n.seq)
    n.min <- min(n.seq)
    max.col <- length(n.seq)
    if (n.min == 1 & is.na(stat(1))) 
        stop("Statistic in stat cannot be computed with n = 1, revise n.seq.")
    if (anim == FALSE) {
        s.dist <- matrix(ncol = 1, nrow = R)
        for (i in 1:R) {
            if(is.expression(parent)) s.dist[i] <- stat(eval(parent))
            if(!is.expression(parent)) s.dist[i] <- stat(sample(parent, size = n.max, replace = FALSE))
        }
        if (!is.null(func)) 
            s.dist <- func(s.dist)
        SE <- round(apply(s.dist,2,sd), digits)
        hist(s.dist, xlab = xlab, main = "", freq = FALSE, 
            breaks = nbreaks, ...)
        if (est.density == TRUE) {
            lines(density(s.dist), col = col.density, lwd = lwd.density)
        }
        if (show.n == TRUE & show.SE == FALSE) {
            legend("topright", legend = paste("n = ", n.max), 
                bty = "n")
        }
        if (show.SE == TRUE & show.n == FALSE) {
            legend("topright", legend = paste("SE = ", SE), bty = "n")
        }
        if (show.SE == TRUE & show.n == TRUE) {
            legend("topright", legend = c(paste("n = ", n.max), 
                paste("SE = ", SE)), bty = "n")
        }
    }
    if (anim == TRUE) {
        if (is.null(col.anim)) 
            clr = NULL
        if (!is.null(col.anim)) {
            if (col.anim == "heat.colors") {
                c1 <- heat.colors(max.col)
                clr <- c1[length(c1):1]
            }
            if (col.anim == "rainbow") 
                clr <- rainbow(max.col)
            if (col.anim == "gray") 
                clr <- gray(1 - (1:max.col/max.col))
            if (col.anim != "rainbow" & col.anim != "gray" & 
                col.anim != "heat.colors") 
                clr <- rep(col.anim, max.col)
        }
        s.dist <- matrix(nrow = R, ncol = max.col)
        for (i in n.min:max.col) {
            if(is.expression(parent)) s.dist[,i] <- stat(eval(parent))
            if(!is.expression(parent)){
            s.dist[,i] <- apply(matrix(replicate(R, sample(parent, 
                i)), i), 2, stat)}
        }
        if (!is.null(func)) 
            s.dist <- func(s.dist)
        SE <- round(apply(s.dist, 2, sd), digits)
        if (est.ylim == TRUE & !is.numeric(ylim)) {
            ylim <- c(0, max(c(max(density(s.dist[, 1])$y), max(density(s.dist[, 
                max.col])$y))))
        }
        brks <- seq(min(s.dist), max(s.dist), length.out = nbreaks)
        for (i in 1:length(n.seq)) {
            dev.hold()
            hist(s.dist[, i], xlab = xlab, ylim = ylim, 
                main = "", breaks = brks, freq = FALSE, col = clr[i],...)
            if (show.n == TRUE & show.SE == FALSE) {
                legend("topright", legend = paste("n = ", n.seq[i]), 
                  bty = "n")
            }
            if (show.SE == TRUE & show.n == FALSE) {
                legend("topright", legend = paste("SE = ", SE[i]), 
                  bty = "n")
            }
            if (show.SE == TRUE & show.n == TRUE) {
                legend("topright", legend = c(paste("n = ", n.seq[i]), 
                  paste("SE = ", SE[i])), bty = "n")
            }
            if (est.density == TRUE) {
                lines(density(s.dist[, i]), col = col.density, 
                  lwd = lwd.density)
            }
            dev.flush()
            Sys.sleep(interval)
        }
    }
}