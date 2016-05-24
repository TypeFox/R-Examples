priceFeature <- function(n, which=c("2clust", "3clust", "5clust",
                            "ellipse", "triangle", "circle", "square",
                            "largesmall"))
{
    which <- match.arg(which)

    circle <- function(n){
        a <- 2*pi*runif(n)
        b <- sqrt(runif(n))
        b * cbind(cos(a),sin(a))
    }

    square <- function(n){
        2 * cbind(runif(n), runif(n)) - 1
    }
        
    
    if(which=="2clust"){
        x <- 2 * circle(round(n/2)) + 3
        x <- rbind(x, 2 * circle(n-nrow(x)) + 7)
    }
    else if(which=="3clust"){
        x <- rbind(1.5 * circle(round(0.4*n)) + 2,
                   1.5 * circle(round(0.4*n)) + 8)
        x <- rbind(x, circle(n-nrow(x)) + 5)
    }
    else if(which=="5clust"){
        x <- rbind(1.5 * circle(round(0.3*n)) + 2,
                   1.5 * circle(round(0.3*n)) + 8)
        x2 <- circle(round(0.12*n)) + 2
        x2[,1] <- x2[,1] + 6
        x3 <- circle(round(0.12*n)) + 2
        x3[,2] <- x2[,2] + 6
        x <- rbind(x,x2,x3)
        x <- rbind(x, circle(n-nrow(x)) + 5)
    }
    else if(which=="ellipse"){
        a <- matrix(c(1,0.5,0.5,1), nrow=2)
        x <- 4 * circle(n) %*% a + 5
    }
    else if(which=="circle"){
        x <- 5 * circle(n) + 5
    }
    else if(which=="square"){
        x <- 4 * square(n) + 5
    }
    else if(which=="triangle"){
        x <- matrix()
        k <- 2.5
        while(nrow(x)<n){
            x <- square(k*n)
            x <- x[x[,1]>x[,2],]
            k <- 1.1*k
        }
        x <- 4*x + 5
    }
    else if(which=="largesmall"){
        x <- circle(round(n/4)) + 8
        x <- rbind(x, 4 * circle(n-nrow(x)) + 4)
    }
    
    x <- x[sample(1:nrow(x)),]    
    class(x) <- "priceFeature"
    colnames(x) <- c("features / performance / quality", "price")
    x
}

plot.priceFeature <- function(x, add=FALSE, ...)
{
    if(add)
        points(x, ...)
    else
        plot.default(x, xlim=c(0,10), ylim=c(0,10), ...)
}
