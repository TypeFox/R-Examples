filtered <-
function(formula, data = NULL, id, process = "bm", timeVar, estimate, subj.id){

mf        <- model.frame(formula = formula, data = data)
resp.comb <- cbind(id, as.matrix(model.extract(mf, "response")))
cov.comb  <- cbind(id, as.matrix(model.matrix(attr(mf, "terms"), data = mf)))
time.comb <- cbind(id, timeVar)   

###############################################
############ BROWNIAN MOTION ##################
###############################################

if(process == "bm"){

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) -1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
tausq.hat   <- estimate[ncol(cov.comb) + 2]

u.save <- w.save <- NULL

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- as.matrix(time.comb[time.comb[, 1] == ii, -1])
ni   <- nrow(X)

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * outer(c(time), c(time), function(x,y) pmin(x, y)) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))


## PREDICTION OF W

for (i in 1:ni){

Fik     <- outer(time[1:i], time[i], function(x,y) pmin(x, y))
Vik.inv <- solve(omegasq.hat * matrix(1, i, i) + 
           sigmasq.hat * outer(time[1:i], time[1:i], function(x,y) pmin(x, y)) + 
           tausq.hat * diag(i))
w.mean  <- sigmasq.hat * t(Fik) %*% Vik.inv %*% (matrix(Y[1:i,]) - X[1:i, ] %*% alpha.hat)
w.var   <- sigmasq.hat * (time[i] - sigmasq.hat * t(Fik) %*% Vik.inv %*% Fik)
w.save  <- rbind(w.save, c(ii, time[i], w.mean, w.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Filtering for the mixed model with Brownian motion"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save

}#bm

##########################################################
############ INTEGRATED BROWNIAN MOTION ##################
##########################################################

if(process == "ibm"){

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) - 1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
tausq.hat   <- estimate[ncol(cov.comb) + 2]

u.save <- b.save <- w.save <- NULL

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- as.matrix(time.comb[time.comb[, 1] == ii, -1])
ni   <- nrow(X)

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * outer(c(time), c(time), function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y))) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var)) 

## PREDICTION OF W

for (i in 1:ni){

Fik     <- outer(time[1:i], time[i], function(x,y) pmin(x, y)^2 * (pmax(x, y) - pmin(x, y) / 3))
Vik.inv <- solve(omegasq.hat * matrix(1, i, i) + 
                 sigmasq.hat * outer(time[1:i], time[1:i], function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))+ 
                 tausq.hat * diag(i))
w.mean  <- 0.5 * sigmasq.hat * t(Fik) %*% Vik.inv %*% (matrix(Y[1:i,]) - X[1:i, ] %*% alpha.hat)
w.var   <- sigmasq.hat * (time[i]^3 / 3 - 0.25 * sigmasq.hat * t(Fik) %*% Vik.inv %*% Fik)
w.save  <- rbind(w.save, c(ii, time[i], w.mean, w.var))

}#i

## PREDICTION OF B

for (i in 1:ni){

Lik     <- matrix(time[1:i]^2, nrow=1)
Vik.inv <- solve(omegasq.hat * matrix(1, i, i) + 
           sigmasq.hat * outer(time[1:i], time[1:i], function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y)))+ 
           tausq.hat * diag(i))

b.mean <- 0.5 * sigmasq.hat * Lik %*% Vik.inv %*% (matrix(Y[1:i,]) - X[1:i, ] %*% alpha.hat)
b.var  <- sigmasq.hat * (time[i] - 0.25 * sigmasq.hat * Lik %*% Vik.inv %*% t(Lik))

b.save <- rbind(b.save, c(ii, time[i], b.mean, b.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- colnames(b.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Filtering for the mixed model with integrated Brownian motion"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save
output$b     <- b.save

}#ibm

#############################################################
############ INTEGRATED ORNSTEIN UHLENBECK ##################
#############################################################

if(process == "iou"){

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) -1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
nu.hat      <- estimate[ncol(cov.comb) + 2]
tausq.hat   <- estimate[ncol(cov.comb) + 3]

u.save <- w.save <- b.save <- c()

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- as.matrix(time.comb[time.comb[, 1] == ii, -1])
ni   <- nrow(X)

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
                outer(c(time), c(time), function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y)))) + 
                tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * matrix(1, 1, ni) %*% Vi.inv %*% t(Ki))

u.save <- rbind(u.save, c(ii, u.mean, u.var))


## PREDICTION OF W

for (i in 1:ni){

Fik     <- outer(time[1:i], time[i], function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y))))
Vik.inv <- solve(omegasq.hat * matrix(1, i, i) + 
                 outer(time[1:i], time[1:i], function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y)))) + 
                 tausq.hat * diag(i))
w.mean <- t(Fik) %*% Vik.inv %*% (matrix(Y[1:i,]) - X[1:i, ] %*% alpha.hat)
w.var  <- sigmasq.hat/(nu.hat^3)*(nu.hat*time[i]+exp(-nu.hat*time[i])-1) - t(Fik) %*% Vik.inv %*% Fik
w.save <- rbind(w.save, c(ii, time[i], w.mean, w.var))

}#i


## PREDICTION OF B

for (i in 1:ni){

Fik     <- outer(time[1:i], time[i], function(x,y) 0.5*sigmasq.hat/(nu.hat^2)*exp(-nu.hat*pmax(x,y))*(exp(nu.hat*pmin(x,y))-1))
Vik.inv <- solve(omegasq.hat * matrix(1, i, i) + 
             outer(time[1:i], time[1:i], function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y)))) + 
             tausq.hat * diag(i))
b.mean <- t(Fik) %*% Vik.inv %*% (matrix(Y[1:i,]) - X[1:i, ] %*% alpha.hat)
b.var  <- 0.5*sigmasq.hat/nu.hat - t(Fik) %*% Vik.inv %*% Fik
b.save <- rbind(b.save, c(ii, time[i], b.mean, b.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- colnames(b.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Filtering for the mixed model with integrated Ornstein-Uhlenbeck process"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save
output$b     <- b.save

}#iou

###########################################
######                              #######
###### STATIONARY GAUSSIAN PROCESS  #######
######  POWERED CORRELATION         #######
###########################################

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){

pow <- as.numeric(strsplit(process, "-")[[1]][3])

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) -1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
phi.hat     <- estimate[ncol(cov.comb) + 2] 
tausq.hat   <- estimate[ncol(cov.comb) + 3]

u.save <- w.save <- NULL

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- as.matrix(time.comb[time.comb[, 1] == ii, -1])
ni   <- nrow(X)

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * outer(c(time), c(time), function(x,y) exp(-abs(x - y)^pow/phi.hat)) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))


## PREDICTION OF W

for (i in 1:ni){

Fik     <- sigmasq.hat * outer(time[1:i], time[i], function(x,y) exp(-abs(x - y)^pow/phi.hat))
Vik.inv <- solve(omegasq.hat * matrix(1, i, i) + 
           sigmasq.hat * outer(time[1:i], time[1:i], function(x,y) exp(-abs(x - y)^pow/phi.hat)) + 
           tausq.hat * diag(i))
w.mean  <- t(Fik) %*% Vik.inv %*% (matrix(Y[1:i,]) - X[1:i, ] %*% alpha.hat)
w.var   <- sigmasq.hat - t(Fik) %*% Vik.inv %*% Fik
w.save  <- rbind(w.save, c(ii, time[i], w.mean, w.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Filtering for the mixed model with Brownian motion"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save
output

}#sgp-powered


###########################################
######                              #######
###### STATIONARY GAUSSIAN PROCESS  #######
######  MATERN CORRELATION          #######
###########################################

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){

kappa <- as.numeric(strsplit(process, "-")[[1]][3])

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) -1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
phi.hat     <- estimate[ncol(cov.comb) + 2] 
tausq.hat   <- estimate[ncol(cov.comb) + 3]

u.save <- w.save <- NULL

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- as.matrix(time.comb[time.comb[, 1] == ii, -1])
ni   <- nrow(X)

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * matern(outer(c(time), c(time), function(x, y) abs(x - y)), phi.hat, kappa = kappa) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))

## PREDICTION OF W

for (i in 1:ni){

Fik     <- sigmasq.hat * matern(outer(time[1:i], time[i], function(x,y) abs(x - y)), phi.hat, kappa = kappa)
Vik.inv <- solve(omegasq.hat * matrix(1, i, i) + 
           sigmasq.hat * matern(outer(time[1:i], time[1:i], function(x, y) abs(x - y)), phi.hat, kappa = kappa) + 
           tausq.hat * diag(i))
w.mean  <- t(Fik) %*% Vik.inv %*% (matrix(Y[1:i,]) - X[1:i, ] %*% alpha.hat)
w.var   <- sigmasq.hat - t(Fik) %*% Vik.inv %*% Fik
w.save  <- rbind(w.save, c(ii, time[i], w.mean, w.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Filtering for the mixed model with stationary Gaussian process - Matern correlation function"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save

}#sgp-matern

output

}
