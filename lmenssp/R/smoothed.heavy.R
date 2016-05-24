smoothed.heavy <-
function(formula, data, id, process, timeVar, estimate, subj.id = NULL){

mf        <- model.frame(formula = formula, data = data)
resp.comb <- cbind(id, as.matrix(model.extract(mf, "response")))
cov.comb  <- cbind(id, as.matrix(model.matrix(attr(mf, "terms"), data = mf)))
time.comb <- cbind(id, timeVar)   

theta <- estimate[, 1]

ncov <- ncol(cov.comb) - 1

alpha <- theta[1 : ncov]
phi   <- theta[(ncov + 1) : (length(theta) - 1)]
dof   <- theta[length(theta)]

u.save <- w.save <- NULL

for (ii in subj.id){

Xi    <- matrix(cov.comb[cov.comb[, 1] == ii, -1], ncol = ncol(cov.comb)-1)
Yi    <- as.matrix(resp.comb[resp.comb[, 1] == ii, -1])
Timei <- time.comb[time.comb[, 1] == ii, -1]
ni   <- nrow(Xi)

Ji <- matrix(1, ni, ni)
Ii <- diag(ni)

if(process == "bm"){
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y)) 
Vi <- phi[1] * Ji + phi[2] * Ri + phi[3] * Ii
}

if(process == "ibm"){
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))   
Vi <- phi[1] * Ji + phi[2] * Ri + phi[3] * Ii
}

if(process == "iou"){
Ri <- outer(Timei, Timei, function(x,y) 
0.5*phi[2]*phi[3]^(-3)*(2*phi[3]*pmin(x,y)+exp(-phi[3]*pmax(x,y))+exp(-phi[3]*pmin(x,y))-1-exp(-phi[3]*abs(pmax(x,y)-pmin(x,y))))
)
Vi <- phi[1] * Ji + Ri + phi[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){
pow <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi[3]))  
Vi <- phi[1] * Ji + phi[2] * Ri + phi[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){
kappa <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- matern(abs(outer(Timei, Timei, "-")), phi[3], kappa = kappa)  
Vi <- phi[1] * Ji + phi[2] * Ri + phi[4] * Ii
}

Vi.inv <- solve(Vi)

deltasq.common <- Yi - Xi %*% alpha
deltasq        <- t(deltasq.common) %*% Vi.inv %*% deltasq.common

eps.hat <- (dof + ni)/(dof + deltasq)

psi.hat <- matrix(0, (ni + 1), (ni + 1))
psi.hat[1, 1] <- phi[1]

if(process == "bm"){
psi.hat[2 : nrow(psi.hat), 2 : nrow(psi.hat)] <- phi[2] * outer(Timei, Timei, function(x,y) pmin(x, y))
}

if(process == "ibm"){
psi.hat[2 : nrow(psi.hat), 2 : nrow(psi.hat)] <- phi[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))
}

if(process == "iou"){
psi.hat[2 : nrow(psi.hat), 2 : nrow(psi.hat)] <- 
outer(Timei, Timei, function(x,y) 
0.5*phi[2]*phi[3]^(-3)*(2*phi[3]*pmin(x,y)+exp(-phi[3]*pmax(x,y))+exp(-phi[3]*pmin(x,y))-1-exp(-phi[3]*abs(pmax(x,y)-pmin(x,y)))))
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){
pow <- as.numeric(strsplit(process, "-")[[1]][3])
psi.hat[2 : nrow(psi.hat), 2 : nrow(psi.hat)] <- phi[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi[3]))
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){
kappa <- as.numeric(strsplit(process, "-")[[1]][3])
psi.hat[2 : nrow(psi.hat), 2 : nrow(psi.hat)] <- phi[2] * matern(abs(outer(Timei, Timei, "-")), phi[3], kappa = kappa)
}

Di <- cbind(matrix(1, ni, 1), diag(ni))
Di.trans <- t(Di)

mean <- as.numeric(psi.hat %*% Di.trans %*% solve(Di %*% psi.hat %*% Di.trans + phi[length(phi)] * diag(ni)) %*% (Yi - Xi %*% matrix(alpha)))
var  <- diag((psi.hat - psi.hat %*% Di.trans %*% solve(Di %*% psi.hat %*% Di.trans + phi[length(phi)] * diag(ni)) %*% Di %*% psi.hat)/
                as.numeric(eps.hat))

u.save <- rbind(u.save, cbind(ii, mean[1], var[1]))
w.save <- rbind(w.save, cbind(ii, Timei, mean[-1], var[-1]))

}##ii 

colnames(u.save) <- c("id", "mean", "variance")
colnames(w.save) <- c("id", "time", "mean", "variance")

output       <- list()
if(process == "bm") output$title <- "Smoothing for the mixed model with t distribution and Brownian motion"
if(process == "ibm") output$title <- "Smoothing for the mixed model with t distribution and integrated Brownian motion"
if(process == "iou") output$title <- "Smoothing for the mixed model with t distribution and integrated Ornstein-Uhlenbeck process"
if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered") 
   output$title <- "Smoothing for the mixed model with t distribution and stationary Gaussian process with powered correlation function"
if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern") 
   output$title <- "Smoothing for the mixed model with t distribution and stationary Gaussian process with Matern correlation function"
output$date  <- date()
output$u     <- u.save
output$w     <- w.save

output

}
