mcmc.norm.hier <- function(data, length = 1000, n.chains = 5){
# require(geoR) 

n.j <- apply(data, 2, function(x) length(x[!is.na(x)]))
ybar.j <- apply(data, 2, function(x) mean(x, na.rm = TRUE))
n <- sum(n.j)
J <- ncol(data)

#--------------- Statistics -----------------#

s.sq.hat <- function(t.j){
sm <- seq(1, J)
for(i in 1 : J){
sm[i] <- sum(sapply(data[,i],function(x)(x - t.j[i])^2),na.rm=TRUE)
}
1/n*sum(sm)
}

tau.sq.hat <- function(t.j, mu) 1/(J-1) * sum((t.j-mu)^2, na.rm = TRUE)
theta.j.hat <- function(tau.sq, s.sq, mu, ybar.j, n.j) ((1/tau.sq) * mu + (n.j/s.sq) * ybar.j)/((1/tau.sq) + (n.j/s.sq))
V.theta.j <- function(tau.sq, s.sq, n.j) 1/((1/tau.sq) + (n.j/s.sq))


#############################################################
#---------------------------Gibbs---------------------------#
#############################################################

dimnames <- list(seq(1 : length), c(paste("theta", seq(1 : J), sep = ""), "mu", "s.sq", "tau.sq"), paste("Chain#", seq(1 : n.chains)))
chains <- array(NA, c(length, J + 3 , n.chains), dimnames = dimnames)
    
    for(m in 1 : n.chains){

#------------------ Start for each chain -------------------#

    theta.j <- matrix(ncol = J, nrow = length)
    mu <- matrix(ncol = 1, nrow = length)
    s.sq <- matrix(ncol = 1, nrow = length)
    tau.sq <- matrix(ncol = 1, nrow = length) 
    
    theta.j[1,] <- apply(data, 2, function(x) sample(x[!is.na(x)], 1, replace = TRUE))
    mu[1] <- mean(theta.j, na.rm = TRUE)
    s.sq1 <- s.sq.hat(theta.j[1,]) 
    s.sq[1] <- rinvchisq(1, n, s.sq1)
    tau.sq1 <- tau.sq.hat(theta.j[1,], mu[1])
    tau.sq[1] <- rinvchisq(1, J - 1, tau.sq1)

#-------------------------- Sampling ------------------------#

    for(i in 2 : length){
        k = i
        for(j in 1 : J){
            theta1 <- theta.j.hat(tau.sq[k - 1], s.sq[k - 1], mu[k - 1],  ybar.j[j], n.j[j])
            vtheta1 <- V.theta.j(tau.sq[k - 1], s.sq[k - 1], n.j[j]) 
            theta.j[k,][j] <- rnorm(1, theta1, sqrt(vtheta1))
            }
    
        mu.hat <- mean(theta.j[i,])
        mu[i] <- rnorm(1, mu.hat, sqrt(tau.sq[i - 1]/(J)))
        
        s.sq1 <- s.sq.hat(theta.j[i,]) 
        s.sq[i] <- rinvchisq(1, n, s.sq1)
        
        tau.sq1 <- tau.sq.hat(theta.j[i,], mu[i])
        tau.sq[i] <- rinvchisq(1, J - 1, tau.sq1)
        } 
    chains[1 : length, 1 : (ncol(data) + 3), m] <- cbind(theta.j, mu, s.sq, tau.sq)
    }

#############################################################
#-----------------------------------------------------------#
#############################################################

class(chains) <- "norm.hier"
chains
}

norm.hier.summary <- function(M, burn.in = 0.5, cred = 0.95, conv.log = TRUE){
    if(class(M) != "norm.hier") stop("Input must be output from function 'mcmc.norm.hier'")
        n <- dim(M)[1]
        m <- M[(n * burn.in) : n, ,]
        vars <- dim(m)[2]
        n.chains <- dim(m)[3]
        
        upper <- 1 - (1 - cred)/2; lower <- 1 - upper 
        res <- matrix(nrow = vars, ncol = 4, dimnames = list(c(paste("theta", seq(1 : (vars - 3)), sep = ""), "mu", "s.sq", "tau.sq"), c(paste(lower * 100, "%", sep = ""), "Median", paste(upper * 100, "%", sep =""), "R.hat")))
            for(i in 1 : vars){
            res[,1][i] <- quantile(m[, i, 1 : n.chains], lower)
            res[,2][i] <- median(m[, i, 1 : n.chains])
            res[,3][i] <- quantile(m[, i, 1 : n.chains], upper)
            }
            if(conv.log == TRUE) M[, (vars - 1) : vars,]  <- log(M[,(vars-1):vars,])
            for(i in 1 : vars){
            res[,4][i] <- R.hat(M[, i, 1 : n.chains])
            }
    res
}


