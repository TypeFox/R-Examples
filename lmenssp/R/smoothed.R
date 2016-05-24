smoothed <-
function(formula, data = NULL, id, process = "bm", timeVar, estimate, subj.id = NULL,  
                     fine = NULL, eq.forec = NULL, uneq.forec = NULL){

mf        <- model.frame(formula = formula, data = data)
resp.comb <- cbind(id, as.matrix(model.extract(mf, "response")))
cov.comb  <- cbind(id, as.matrix(model.matrix(attr(mf, "terms"), data = mf)))
time.comb <- cbind(id, timeVar)   

if(length(fine) == 0 & length(eq.forec) == 0 & length(uneq.forec) == 0){
time.comb.sub <- matrix(time.comb[time.comb[, 1] %in% subj.id, ], ncol = 2)
timeSmooth.all <- tapply(time.comb.sub[, 2], time.comb.sub[, 1], function(x) x)
} else if(length(fine) > 0 & length(eq.forec) == 0 & length(uneq.forec) == 0){
time.comb.sub <- matrix(time.comb[time.comb[, 1] %in% subj.id, ], ncol = 2)
timeSmooth.all <- tapply(time.comb.sub[, 2], time.comb.sub[, 1], function(x) seq(min(x), round(max(x), (nchar(fine) - 2)), by = fine))
} else if(length(fine) == 0 & length(eq.forec) > 0 & length(uneq.forec) == 0){
time.comb.sub <- matrix(time.comb[time.comb[, 1] %in% subj.id, ], ncol = 2)
timeSmooth.all <- tapply(time.comb.sub[, 2], time.comb.sub[, 1], function(x) max(x) + cumsum(rep(eq.forec[1], eq.forec[2])))
} else {
subj.id        <- unique(uneq.forec[, 1])
nobs           <- as.numeric(table(uneq.forec[, 1]))
time.comb.sub  <- time.comb[time.comb[, 1] %in% subj.id, ]
max.time       <- as.numeric(unlist(tapply(time.comb.sub[, 2], time.comb.sub[, 1], function(x) max(x))))
max.time.ext   <- rep(max.time, nobs)
uneq.forec[, 2] <- max.time.ext + uneq.forec[, 2]   
timeSmooth.all <- tapply(uneq.forec[, 2], uneq.forec[, 1], function(x) x)
}


###############################################
############ BROWNIAN MOTION ##################
###############################################

if(process == "bm"){

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) - 1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
tausq.hat   <- estimate[ncol(cov.comb) + 2]

u.save <- w.save <- NULL

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- time.comb[time.comb[, 1] == ii, -1]
ni   <- nrow(X)

timeSmooth <- unlist(timeSmooth.all[which((subj.id==ii) == TRUE)])

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * outer(c(time), c(time), function(x,y) pmin(x, y)) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))

## SMOOTHING W

for(i in 1 : length(timeSmooth)){

Fi     <- outer(timeSmooth[i], time, function(x,y) pmin(x, y))
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
                sigmasq.hat * outer(c(time), c(time), function(x,y) pmin(x, y)) + 
                tausq.hat * diag(ni))
w.mean <- sigmasq.hat * Fi %*% Vi.inv %*% (Y - X %*% alpha.hat)
w.var  <- sigmasq.hat * (timeSmooth[i] - sigmasq.hat * Fi %*% Vi.inv %*% t(Fi))
w.save <- rbind(w.save, c(ii, timeSmooth[i], w.mean, w.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Smoothing for the mixed model with Brownian motion"
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

u.save <- w.save <- b.save <- NULL

for(ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- time.comb[time.comb[, 1] == ii, -1]
ni   <- nrow(X)

timeSmooth <- unlist(timeSmooth.all[which((subj.id==ii) == TRUE)])

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * outer(c(time), c(time), function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y))) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))

## SMOOTHING W

for(i in 1 : length(timeSmooth)){

Fi     <- outer(timeSmooth[i], time, function(x,y) pmin(x, y)^2 * (pmax(x, y) - pmin(x, y) / 3)) 
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
                sigmasq.hat * outer(c(time), c(time), function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y))) + 
                tausq.hat * diag(ni))
w.mean <- 0.5 * sigmasq.hat * Fi %*% Vi.inv %*% (Y - X %*% alpha.hat) 
w.var  <- sigmasq.hat * (timeSmooth[i]^3/3 - 0.25 * sigmasq.hat * Fi %*% Vi.inv %*% t(Fi))
w.save <- rbind(w.save, c(ii, timeSmooth[i], w.mean, w.var))

}#i

## SMOOTHING B

for(i in 1 : length(timeSmooth)){

L1     <- 0.5 * sigmasq.hat * outer(timeSmooth[i], time, function(x,y) pmin(x, y)^2)
L2     <- sigmasq.hat * outer(timeSmooth[i], time, function(x,y) x*y - 0.5*x^2)
Li     <- matrix(c(L1[1, timeSmooth[i] >= time], L2[1, timeSmooth[i] < time]), nrow = 1)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
                sigmasq.hat * outer(c(time), c(time), function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - 1/3 * pmin(x, y))) + 
                tausq.hat * diag(ni))
b.mean <- Li %*% Vi.inv %*% (Y - X %*% alpha.hat)
b.var  <- sigmasq.hat * timeSmooth[i] - Li %*% Vi.inv %*% t(Li)
b.save <- rbind(b.save, c(ii, timeSmooth[i], b.mean, b.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- colnames(b.save) <- c("id", "time", "mean", "variance") 

output       <- list()
output$title <- "Smoothing for the mixed model with integrated Brownian motion"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save
output$b     <- b.save
output

}

#############################################################
############ INTEGRATED ORNSTEIN UHLENBECK ##################
#############################################################

if(process == "iou"){

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) - 1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
nu.hat      <- estimate[ncol(cov.comb) + 2]
tausq.hat   <- estimate[ncol(cov.comb) + 3]

u.save <- w.save <- b.save <- NULL

for(ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- time.comb[time.comb[, 1] == ii, -1]
ni   <- nrow(X)

timeSmooth <- unlist(timeSmooth.all[which((subj.id==ii) == TRUE)])

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
                outer(c(time), c(time), function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y)))) + 
                tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * matrix(1, 1, ni) %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))

## SMOOTHING W

for(i in 1 : length(timeSmooth)){

F <- outer(timeSmooth[i], time, function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y))))
V.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
               outer(time, time, function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y)))) + 
               tausq.hat * diag(ni)) 
w.mean <- F %*% V.inv %*% (Y - X %*% alpha.hat)
w.var  <- sigmasq.hat/(nu.hat^3)*(nu.hat*timeSmooth[i]+exp(-nu.hat*timeSmooth[i])-1) - round(F %*% V.inv %*% t(F), 20)
w.save <- rbind(w.save, c(ii, timeSmooth[i], w.mean, w.var))

}#i

## SMOOTHING B

b.save <- NULL

for(i in 1 : length(timeSmooth)){

F1    <- outer(timeSmooth[i], time, function(x,y) 0.5*sigmasq.hat/(nu.hat^2)*exp(-nu.hat*x)*(exp(nu.hat*y)-1))
F2    <- outer(timeSmooth[i], time, function(x,y) 0.5*sigmasq.hat/(nu.hat^2)*(2-exp(-nu.hat*x)-exp(nu.hat*(x-y))) )
F     <- matrix(c(F1[1, timeSmooth[i] >= time], F2[1, timeSmooth[i] < time]), nrow = 1)
V.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
               outer(time, time, function(x,y) 0.5*sigmasq.hat/(nu.hat^3)*(2*nu.hat*pmin(x,y)+exp(-nu.hat*x)+exp(-nu.hat*y)-1-exp(-nu.hat*abs(x-y)))) + 
               tausq.hat * diag(ni))
b.mean <- F %*% V.inv %*% (Y - X %*% alpha.hat)
b.var  <- 0.5*sigmasq.hat/nu.hat - F %*% V.inv %*% t(F)
b.save <- rbind(b.save, c(ii, timeSmooth[i], b.mean, b.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- colnames(b.save) <- c("id", "time", "mean", "variance") 

output       <- list()
output$title <- "Smoothing for the mixed model with integrated Ornstein-Uhlenbeck process"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save
output$b     <- b.save

}#iou

########################################
#####                             ###### 
##### STATIONARY GAUSSIAN PROCESS ######
#####   POWERED CORRELATION       ######
########################################
 
if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){

pow <- as.numeric(strsplit(process, "-")[[1]][3])

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) - 1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
phi.hat     <- estimate[ncol(cov.comb) + 2]
tausq.hat   <- estimate[ncol(cov.comb) + 3]

u.save <- w.save <- NULL

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- time.comb[time.comb[, 1] == ii, -1]
ni   <- nrow(X)

timeSmooth <- unlist(timeSmooth.all[which((subj.id==ii) == TRUE)])

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * outer(c(time), c(time), function(x,y) exp(-abs(x - y)^pow/phi.hat)) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))

## SMOOTHING W

for(i in 1 : length(timeSmooth)){

Fi     <- sigmasq.hat * outer(timeSmooth[i], time, function(x,y) exp(-abs(x - y)^pow/phi.hat))
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
                sigmasq.hat * outer(c(time), c(time), function(x,y) exp(-abs(x - y)^pow/phi.hat)) + 
                tausq.hat * diag(ni))
w.mean <- Fi %*% Vi.inv %*% (Y - X %*% alpha.hat)
w.var  <- sigmasq.hat - Fi %*% Vi.inv %*% t(Fi)
w.save <- rbind(w.save, c(ii, timeSmooth[i], w.mean, w.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Smoothing for the mixed model with stationary Gaussian process - powered correlation function"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save
output

}#sgp-powered

########################################
#####                             ###### 
##### STATIONARY GAUSSIAN PROCESS ######
#####   MATERN CORRELATION        ######
########################################
 
if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){

kappa <- as.numeric(strsplit(process, "-")[[1]][3])

alpha.hat   <- matrix(estimate[1 : (ncol(cov.comb) - 1)])
omegasq.hat <- estimate[ncol(cov.comb)]
sigmasq.hat <- estimate[ncol(cov.comb) + 1]
phi.hat     <- estimate[ncol(cov.comb) + 2]
tausq.hat   <- estimate[ncol(cov.comb) + 3]

u.save <- w.save <- NULL

for (ii in subj.id){

X    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Y    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
time <- time.comb[time.comb[, 1] == ii, -1]
ni   <- nrow(X)

timeSmooth <- unlist(timeSmooth.all[which((subj.id==ii) == TRUE)])

## PREDICTION OF U

Ki     <- matrix(1, 1, ni)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
          sigmasq.hat * matern(u = outer(c(time), c(time), function(x,y) abs(x-y)), phi = phi.hat, kappa = kappa) + 
          tausq.hat * diag(ni))
u.mean <- omegasq.hat * Ki %*% Vi.inv %*% (Y - X %*% alpha.hat)
u.var  <- omegasq.hat * (1 - omegasq.hat * Ki %*% Vi.inv %*% t(Ki))
u.save <- rbind(u.save, c(ii, u.mean, u.var))

## SMOOTHING W

for(i in 1 : length(timeSmooth)){

Fi     <- sigmasq.hat * matern(u = outer(timeSmooth[i], time, function(x,y) abs(x-y)), phi = phi.hat, kappa = kappa)
Vi.inv <- solve(omegasq.hat * matrix(1, ni, ni) + 
                sigmasq.hat * matern(u = outer(c(time), c(time), function(x,y) abs(x-y)), phi = phi.hat, kappa= kappa) + 
                tausq.hat * diag(ni))
w.mean <- Fi %*% Vi.inv %*% (Y - X %*% alpha.hat)
w.var  <- sigmasq.hat - Fi %*% Vi.inv %*% t(Fi)
w.save <- rbind(w.save, c(ii, timeSmooth[i], w.mean, w.var))

}#i
}#ii

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- c("id", "time", "mean", "variance")

output       <- list()
output$title <- "Smoothing for the mixed model with stationary Gaussian process - Matern correlation function"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save

}#sgp-matern

output

}
