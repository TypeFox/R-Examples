anm.mc.norm <- function(start = -4, mu = 0, sigma = 1, length = 2000, sim = "M", jump.kernel = 0.2, xlim = c(-4, 4), ylim = c(0, 0.4), interval = 0.01, show.leg = TRUE, ...)
{
x<-seq(-4, 4, .05)
dev.new(xpos = -750)
curve(dnorm(x, mu, sigma), xlim = xlim, xlab = "x", ylab = "f(x)", main = "Underlying distribution")
dev.new()

if(sim == "M") method = "Metropolis"; if(sim == "MH") method = "Metropolis-Hastings"

q.x <- matrix(ncol = 1, nrow = length)
q.x[1] <- start

if(sim == "M")
{
    for( i in 2 : length)
        {
        new <- rnorm(1, q.x[i - 1], sigma * (jump.kernel^2))
        dold <- dnorm(q.x[i - 1], mu, sigma)
        dnew <- dnorm(new, mu, sigma)
        r <- dnew/dold
        p <- min(r, 1)
        s <- sample(c(1, 2), 1, prob = c(p, 1 - p))
        if(s == 1)q.x[i] = new
        if(s == 2)q.x[i] = q.x[i - 1]
        }
}

if(sim == "MH")
{        

    for( i in 2 : length)
        {
        new <- rnorm(1, q.x[i - 1], sigma * (jump.kernel^2))
        dold <- dnorm(q.x[i - 1], mu, sigma)
        djump.old <- dnorm(q.x[i - 1], q.x[i - 1], sigma)
        dnew <- dnorm(new, mu, sigma)
        djump.new <- dnorm(new, new, sigma)
        r <- (dnew/djump.new)/(dold/djump.old)
        p <- min(r, 1)
        s <- sample(c(1, 2), 1, prob = c(p, 1 - p))
        if(s == 1)q.x[i] = new
        if(s == 2)q.x[i] = q.x[i - 1]
        }
}        
        
dens <- dnorm(q.x, mu, sigma)
    
    for( i in 1 : length)
        {
        dev.hold()
        plot(q.x[seq(1 : i)], dens[seq(1 : i)], xlim = xlim, ylim = ylim, xlab = "x", ylab = "f(x)", type = "l", col = "gray", main = paste("MCMC walk, method = ", method))
        points(q.x[seq(1 : i)], dens[seq(1 : i)], pch = 19, cex = .5)
        points(q.x[i], dens[i], pch = 19, cex = 1.5)
        legend("topleft", bty = "n", legend = c("Length = ", i))
        dev.flush()
        Sys.sleep(interval)
        }
}




anm.mc.bvn <- function(start = c(-4, -4), mu = c(0, 0), sigma = matrix(2, 2, data = c(1, 0, 0, 1)), length = 1000, sim = "M", jump.kernel = 0.2, xlim = c(-4, 4), ylim = c(-4, 4), interval = 0.01, show.leg = TRUE, cex.leg = 1,...)
{


x<-seq(-4, 4, .05)
g<-expand.grid(x, x)
p<-dmvnorm(g, mu, sigma)
dev.new(xpos = -750)
old.par <- par(no.readonly = TRUE)
par(mar=c(5,4.5,3,2))

plot(g, type = "p", col = gray(1-(p/max(p))), xlab = expression(italic(X)[1]), ylab = expression(italic(X)[2]), main = "Underlying BVN distribution",...)
dev.new()
par(mar=c(5,4.5,3,2))

if(sim == "M") method = "Metropolis"; if(sim == "MH") method = "Metropolis-Hastings"; if(sim == "G") method = "Gibbs"

q.x <- matrix(ncol = 2, nrow = length)
q.x[1,] <- start

if(sim == "G")
{
x <- y <- matrix(ncol = 1, nrow = length)
x[1] <- start[1]
y[1] <- start[2] 
rho <- sigma[3]
    for( i in 2 : length)
        {
        x[i] <- rnorm(1, mu[1] + rho * y[i - 1], sd = 1 - rho^2)
        y[i] <- rnorm(1, mu[2] + rho * x[i], sd = 1 - rho^2)
        }
        q.x <- cbind(x, y)
}


if(sim == "M")
{
    for( i in 2 : length)
        {
        x <- rmvnorm(1, q.x[i - 1,], sigma * (jump.kernel^2))
        dold <- dmvnorm(q.x[i - 1,], mu, sigma)
        dnew <- dmvnorm(x, mu, sigma)
        r <- dnew/dold
        p <- min(1, r)
        s <- sample(c(1, 2), 1, prob = c(p, 1 - p))
        if(s == 1)q.x[i,] = x
        if(s == 2)q.x[i,] = q.x[i - 1,]
        }
}


if(sim == "MH")
{        
    for( i in 2 : length)
        {
        x <- rmvnorm(1, q.x[i - 1,], sigma * (jump.kernel^2))
        dold <- dmvnorm(q.x[i - 1,], mu, sigma)
        djump.old <- dmvnorm(q.x[i - 1,], q.x[i - 1,], sigma)
        dnew <- dmvnorm(x, mu, sigma)
        djump.new <- dmvnorm(x, x, sigma)
        r <- (dnew/djump.new)/(dold/djump.old)
        p <- min(r, 1)
        s <- sample(c(1, 2), 1, prob = c(p, 1 - p))
        if(s == 1)q.x[i,] = x
        if(s == 2)q.x[i,] = q.x[i - 1,]
        }
}        
        
    for( i in 1 : length)
        {
        dev.hold()
        plot(q.x[seq(1 : i),], xlim = xlim, ylim = ylim, xlab = expression(italic(X)[1]), ylab = expression(italic(X)[2]), type = "l", col = "gray", main = paste("MCMC walk, method = ", method),...)
        points(q.x[,1][seq(1 : i)], q.x[,2][seq(1 : i)], pch = 19, cex = .5)
        points(q.x[i,][1], q.x[i,][2], pch = 19, cex = 1.5)
        density <- dmvnorm(q.x[i,], mu, sigma)
        #legend("topright", bty = "n", inset = 0.1, legend = c("f(x1,x2) = ", density))
        if(show.leg) legend("topleft", bty = "n", legend = c("Length = ", i), cex = cex.leg)
        dev.flush()
        Sys.sleep(interval)
        }
on.exit(par(old.par))        
}
                                                                                                                                            