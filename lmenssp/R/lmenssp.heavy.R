lmenssp.heavy <- function(formula, data, id, timeVar, init.em = NULL, maxiter.em = 100, tol.em = 0.001, process, silent = TRUE, 
                          dof.est = c(0.1, 10, 0.0001), tol.cd = 0.001,
                          tol.lmenssp = 10^-5, init.lmenssp = NULL, maxiter.lmenssp = 100, silent.lmenssp = TRUE){

#library(geoR)
#library(nlme)
#library(mvtnorm)
#library(MASS)

init <- init.em

if(length(init) == 0){
fit.init <- lmenssp(formula = formula, data = data, id = id, timeVar = timeVar, 
               process = process)
init <- c(fit.init$est[, 1], 5)
} 
if(length(init) == 1){
fit.init <- lmenssp(formula = formula, data = data, id = id, timeVar = timeVar, 
               process = process)
init <- c(fit.init$est[, 1], init)
} 


param.new <- init

param.old <- param.new * 20

iter.em   <- 1
Niter.em  <- maxiter.em
tol.em    <- tol.em

id.em <- id
timeVar.em <- timeVar

nobs.em   <- as.numeric(table(id.em))
idlist.em <- unique(id.em)  

formula.em <- formula

mf.em <- model.frame(formula = formula.em, data = data)
y.em  <- as.matrix(model.extract(mf.em, "response"))
colnames(y.em) <- "y"
x.em  <- as.matrix(model.matrix(attr(mf.em, "terms"), data = mf.em))
colnames(x.em)[1] <- gsub("[[:punct:]]", "", colnames(x.em)[1])

nsubj.em  <- length(unique(id.em)) # number of subjects

Time.em  <- tapply(timeVar.em, id.em, function(x) x)
data2.em <- data.frame(cbind(id.em, x.em))
DM.em    <- split(data2.em[, -1], data2.em$id.em)
YM.em    <- tapply(y.em, id.em, function(x) x)
nobs.em  <- as.numeric(tapply(id.em, id.em, function(x) length(x)))

while(sqrt((param.old - param.new) %*% (param.old - param.new)) > tol.em & iter.em <= Niter.em){

param.old <- param.new

alpha <- as.matrix(param.old[1 : ncol(x.em)])
phi   <- param.old[(ncol(x.em) + 1) : (length(param.new) - 1)]
dof   <- param.old[length(param.new)]


### mean of epsilon|data

eps.hat <- c()

for(i in 1 : nsubj.em){

Timei <- Time.em[[i]]
ni    <- nobs.em[i]
Yi    <- as.matrix(YM.em[[i]], nrow = ni)
Xi    <- as.matrix(DM.em[[i]], nrow = ni)

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

eps.hat <- c(eps.hat, (dof + ni)/(dof + deltasq))

}#for(i in idlist.em)

eps.hat.long <- rep(eps.hat, nobs.em)


## maximising wrt theta excluding dof

data2 <- as.data.frame(cbind(y.em * sqrt(eps.hat.long), x.em * sqrt(eps.hat.long), timeVar.em))

formula2 <- reformulate(termlabels = names(data2)[2:(ncol(data2)-1)], response = names(data2)[1], intercept = FALSE)

if(silent.lmenssp == FALSE){

print('---------------------------')
print('------ lmenssp steps ------')
print('---------------------------')

}

fit <- lmenssp(formula = formula2, data = data2, id = id.em, timeVar = timeVar.em, 
               process = process, tol = tol.lmenssp, init = init.lmenssp, maxiter = maxiter.lmenssp, 
               silent = silent.lmenssp)
theta.hat <- fit$est[, 1]

## maximising wrt dof

dof.hat <- function(nu){

alpha.hat <- as.matrix(theta.hat[1 : ncol(x.em)])
phi.hat   <- theta.hat[(ncol(x.em) + 1) : length(theta.hat)]

sum1 <- 0

for (i in 1 : nsubj.em){

Timei <- Time.em[[i]]
ni    <- nobs.em[i]
Yi    <- as.matrix(YM.em[[i]], nrow = ni)
Xi    <- as.matrix(DM.em[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ii <- diag(ni)

if(process == "bm"){
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y)) 
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[3] * Ii
}

if(process == "ibm"){
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))   
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[3] * Ii
}

if(process == "iou"){
Ri <- outer(Timei, Timei, function(x,y) 
0.5*phi.hat[2]*phi.hat[3]^(-3)*(2*phi.hat[3]*pmin(x,y)+exp(-phi.hat[3]*pmax(x,y))+exp(-phi.hat[3]*pmin(x,y))-1-exp(-phi.hat[3]*abs(pmax(x,y)-pmin(x,y))))
)
Vi <- phi.hat[1] * Ji + Ri + phi.hat[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){
pow <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi.hat[3]))  
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){
kappa <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- matern(abs(outer(Timei, Timei, "-")), phi.hat[3], kappa = kappa)  
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[4] * Ii
}

Vi.inv <- solve(Vi)

deltasq.common <- Yi - Xi %*% alpha.hat
deltasq        <- t(deltasq.common) %*% Vi.inv %*% deltasq.common

sum1 <- sum1 + log(gamma((nu + ni) * 0.5)) - log(gamma(nu * 0.5)) + 0.5 * nu * log(nu) - 0.5 * (nu + ni) * log(nu + deltasq) 

}#for

-sum1

}#dof.hat

nu.hat <- optimize(dof.hat, c(dof.est[1], dof.est[2]), tol = dof.est[3])$minimum

param.new <- c(theta.hat, nu.hat)
names(param.new)[length(param.new)] <- "dof"

if(silent == FALSE){

print('--------')
print("EM steps")
print('--------')

cat("iteration = ", iter.em, "\n")
cat("alpha.old = ", param.old[1 : ncol(x.em)], "\n")
cat("alpha.new = ", param.new[1 : ncol(x.em)], "\n")
cat("varpar.old = ", param.old[(ncol(x.em) + 1) : (length(param.old) - 1)], "\n")
cat("varpar.new = ", param.new[(ncol(x.em) + 1) : (length(param.new) - 1)], "\n")
cat("dof.old = ", param.old[length(param.old)], "\n")
cat("dof.new = ", param.new[length(param.new)], "\n")
cat("sqrt.diff = ", sqrt((param.old - param.new) %*% (param.old - param.new)), "\n")
print("-----------------------------")

}

iter.em <- iter.em + 1

}#while


##### standard errors for alpha

phi.hat <- param.new[(ncol(x.em) + 1) : (length(param.new) - 1)]
dof.hat <- param.new[length(param.new)]

hessian <- 0

for(i in 1 : nsubj.em){

Timei <- Time.em[[i]]
ni    <- nobs.em[i]
Xi    <- as.matrix(DM.em[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ii <- diag(ni)

if(process == "bm"){
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y)) 
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[3] * Ii
}

if(process == "ibm"){
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))   
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[3] * Ii
}

if(process == "iou"){
Ri <- outer(Timei, Timei, function(x,y) 
0.5*phi.hat[2]*phi.hat[3]^(-3)*(2*phi.hat[3]*pmin(x,y)+exp(-phi.hat[3]*pmax(x,y))+exp(-phi.hat[3]*pmin(x,y))-1-exp(-phi.hat[3]*abs(pmax(x,y)-pmin(x,y))))
)
Vi <- phi.hat[1] * Ji + Ri + phi.hat[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){
pow <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi.hat[3]))  
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){
kappa <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- matern(abs(outer(Timei, Timei, "-")), phi.hat[3], kappa = kappa)  
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[4] * Ii
}

Vi.inv <- solve(Vi)

hessian <- hessian + (dof.hat + ni)/(dof.hat + ni + 2) * t(Xi) %*% Vi.inv %*% Xi

}#for(i in idlist.em)


#### central difference approximation

theta <- param.new

h <- tol.cd

hessian.full <- 0

for(i in 1 : nsubj.em){

Timei <- Time.em[[i]]
ni    <- nobs.em[i]
Yi    <- as.matrix(YM.em[[i]], nrow = ni)
Xi    <- as.matrix(DM.em[[i]], nrow = ni)

hessian.i <- matrix(0, length(theta), length(theta))

for(j in 1 : length(theta)){
for(k in 1 : length(theta)){

if(j == k){

theta1 <- theta2 <- theta3 <- theta
theta1[j] <- theta1[j] + 2*h
theta2    <- theta
theta3[j] <- theta3[j] - 2*h

alpha1 <- theta1[1 : ncol(x.em)]; phi1 <- theta1[(ncol(x.em) + 1) : (length(theta1) - 1)]; dof1 <- theta1[length(theta1)]
alpha2 <- theta2[1 : ncol(x.em)]; phi2 <- theta2[(ncol(x.em) + 1) : (length(theta2) - 1)]; dof2 <- theta2[length(theta2)]
alpha3 <- theta3[1 : ncol(x.em)]; phi3 <- theta3[(ncol(x.em) + 1) : (length(theta3) - 1)]; dof3 <- theta3[length(theta3)]

if(process == "bm"){

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * outer(Timei, Timei, function(x,y) pmin(x, y)) + 
        phi1[3] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * outer(Timei, Timei, function(x,y) pmin(x, y)) + 
        phi2[3] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * outer(Timei, Timei, function(x,y) pmin(x, y)) + 
        phi3[3] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

}#bm

if(process == "ibm"){

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3)) + 
        phi1[3] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3)) + 
        phi2[3] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3)) + 
        phi3[3] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

}#ibm

if(process == "iou"){

V1   <- phi1[1] * matrix(1, ni, ni) +  
        outer(Timei, Timei, function(x,y) 
          0.5*phi1[2]*phi1[3]^(-3)*(2*phi1[3]*pmin(x,y)+exp(-phi1[3]*pmax(x,y))+exp(-phi1[3]*pmin(x,y))-1-exp(-phi1[3]*abs(pmax(x,y)-pmin(x,y))))) + 
        phi1[3] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        outer(Timei, Timei, function(x,y) 
         0.5*phi2[2]*phi2[3]^(-3)*(2*phi2[3]*pmin(x,y)+exp(-phi2[3]*pmax(x,y))+exp(-phi2[3]*pmin(x,y))-1-exp(-phi2[3]*abs(pmax(x,y)-pmin(x,y))))) + 
        phi2[3] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        outer(Timei, Timei, function(x,y) 
         0.5*phi3[2]*phi3[3]^(-3)*(2*phi3[3]*pmin(x,y)+exp(-phi3[3]*pmax(x,y))+exp(-phi3[3]*pmin(x,y))-1-exp(-phi3[3]*abs(pmax(x,y)-pmin(x,y))))) + 
        phi3[3] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

}#iou

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){

pow <- as.numeric(strsplit(process, "-")[[1]][3])

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi1[3])) + 
        phi1[4] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi2[3])) + 
        phi2[4] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi3[3])) + 
        phi3[4] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

}#powered

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){

kappa <- as.numeric(strsplit(process, "-")[[1]][3])

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * matern(abs(outer(Timei, Timei, "-")), phi1[3], kappa = kappa) + 
        phi1[4] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * matern(abs(outer(Timei, Timei, "-")), phi2[3], kappa = kappa) + 
        phi2[4] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * matern(abs(outer(Timei, Timei, "-")), phi3[3], kappa = kappa) + 
        phi3[4] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

}#matern

hessian.i[j, k] <- (exp1 - 2*exp2 + exp3)/(4*h^2)

}#if(j == k)

if(j != k){

theta1 <- theta2 <- theta3 <- theta4 <- theta
theta1[j] <- theta1[j] + h; theta1[k] <- theta1[k] + h
theta2[j] <- theta2[j] + h; theta2[k] <- theta2[k] - h
theta3[j] <- theta3[j] - h; theta3[k] <- theta3[k] + h
theta4[j] <- theta4[j] - h; theta4[k] <- theta4[k] - h

alpha1 <- theta1[1 : ncol(x.em)]; phi1 <- theta1[(ncol(x.em) + 1) : (length(theta1) - 1)]; dof1 <- theta1[length(theta1)]
alpha2 <- theta2[1 : ncol(x.em)]; phi2 <- theta2[(ncol(x.em) + 1) : (length(theta2) - 1)]; dof2 <- theta2[length(theta2)]
alpha3 <- theta3[1 : ncol(x.em)]; phi3 <- theta3[(ncol(x.em) + 1) : (length(theta3) - 1)]; dof3 <- theta3[length(theta3)]
alpha4 <- theta4[1 : ncol(x.em)]; phi4 <- theta4[(ncol(x.em) + 1) : (length(theta4) - 1)]; dof4 <- theta4[length(theta4)]

if(process == "bm"){

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * outer(Timei, Timei, function(x,y) pmin(x, y)) + 
        phi1[3] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * outer(Timei, Timei, function(x,y) pmin(x, y)) + 
        phi2[3] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * outer(Timei, Timei, function(x,y) pmin(x, y)) + 
        phi3[3] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

V4   <- phi4[1] * matrix(1, ni, ni) + 
        phi4[2] * outer(Timei, Timei, function(x,y) pmin(x, y)) + 
        phi4[3] * diag(ni)
exp4 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha4), sigma = V4, df = dof4, log = TRUE)

}#bm

if(process == "ibm"){

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3)) + 
        phi1[3] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3)) + 
        phi2[3] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3)) + 
        phi3[3] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

V4   <- phi4[1] * matrix(1, ni, ni) + 
        phi4[2] * outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3)) + 
        phi4[3] * diag(ni)
exp4 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha4), sigma = V4, df = dof4, log = TRUE)

}#ibm

if(process == "iou"){

V1   <- phi1[1] * matrix(1, ni, ni) +  
        outer(Timei, Timei, function(x,y) 
         0.5*phi1[2]*phi1[3]^(-3)*(2*phi1[3]*pmin(x,y)+exp(-phi1[3]*pmax(x,y))+exp(-phi1[3]*pmin(x,y))-1-exp(-phi1[3]*abs(pmax(x,y)-pmin(x,y))))) + 
        phi1[4] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        outer(Timei, Timei, function(x,y) 
         0.5*phi2[2]*phi2[3]^(-3)*(2*phi2[3]*pmin(x,y)+exp(-phi2[3]*pmax(x,y))+exp(-phi2[3]*pmin(x,y))-1-exp(-phi2[3]*abs(pmax(x,y)-pmin(x,y))))) + 
        phi2[4] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        outer(Timei, Timei, function(x,y) 
         0.5*phi3[2]*phi3[3]^(-3)*(2*phi3[3]*pmin(x,y)+exp(-phi3[3]*pmax(x,y))+exp(-phi3[3]*pmin(x,y))-1-exp(-phi3[3]*abs(pmax(x,y)-pmin(x,y))))) + 
        phi3[4] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

V4   <- phi4[1] * matrix(1, ni, ni) + 
        outer(Timei, Timei, function(x,y) 
         0.5*phi4[2]*phi4[3]^(-3)*(2*phi4[3]*pmin(x,y)+exp(-phi4[3]*pmax(x,y))+exp(-phi4[3]*pmin(x,y))-1-exp(-phi4[3]*abs(pmax(x,y)-pmin(x,y))))) + 
        phi4[4] * diag(ni)
exp4 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha4), sigma = V4, df = dof4, log = TRUE)

}#iou

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){

pow <- as.numeric(strsplit(process, "-")[[1]][3])

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi1[3])) + 
        phi1[4] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi2[3])) + 
        phi2[4] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi3[3])) + 
        phi3[4] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

V4   <- phi4[1] * matrix(1, ni, ni) + 
        phi4[2] * outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi4[3])) + 
        phi4[4] * diag(ni)
exp4 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha4), sigma = V4, df = dof4, log = TRUE)

}#powered

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){

kappa <- as.numeric(strsplit(process, "-")[[1]][3])

V1   <- phi1[1] * matrix(1, ni, ni) +  
        phi1[2] * matern(abs(outer(Timei, Timei, "-")), phi1[3], kappa = kappa) + 
        phi1[4] * diag(ni)
exp1 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha1), sigma = V1, df = dof1, log = TRUE)

V2   <- phi2[1] * matrix(1, ni, ni) +  
        phi2[2] * matern(abs(outer(Timei, Timei, "-")), phi2[3], kappa = kappa) + 
        phi2[4] * diag(ni)
exp2 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha2), sigma = V2, df = dof2, log = TRUE)

V3   <- phi3[1] * matrix(1, ni, ni) + 
        phi3[2] * matern(abs(outer(Timei, Timei, "-")), phi3[3], kappa = kappa) + 
        phi3[4] * diag(ni)
exp3 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha3), sigma = V3, df = dof3, log = TRUE)

V4   <- phi4[1] * matrix(1, ni, ni) + 
        phi4[2] * matern(abs(outer(Timei, Timei, "-")), phi4[3], kappa = kappa) + 
        phi4[4] * diag(ni)
exp4 <- dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha4), sigma = V4, df = dof4, log = TRUE)

}#kappa

hessian.i[j, k] <- (exp1 - exp2 - exp3 + exp4)/(4*h^2)

}#if(j != k)

}#for(k in 1 : length(theta)
}#for(j in 1 : length(theta)

hessian.full <- hessian.full + hessian.i

#print(i)

}#for(i in idlist)


#### maximised log-likelihood

alpha.hat <- as.matrix(param.new[1 : ncol(x.em)])
phi.hat   <- param.new[(ncol(x.em) + 1) : (length(param.new) - 1)]
nu.hat    <- param.new[length(param.new)]

sum.llik <- 0

for(i in 1 : nsubj.em){

Timei <- Time.em[[i]]
ni    <- nobs.em[i]
Yi    <- as.matrix(YM.em[[i]], nrow = ni)
Xi    <- as.matrix(DM.em[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ii <- diag(ni)

if(process == "bm"){
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y)) 
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[3] * Ii
}

if(process == "ibm"){
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))   
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[3] * Ii
}

if(process == "iou"){
Ri <- outer(Timei, Timei, function(x,y) 
0.5*phi.hat[2]*phi.hat[3]^(-3)*(2*phi.hat[3]*pmin(x,y)+exp(-phi.hat[3]*pmax(x,y))+exp(-phi.hat[3]*pmin(x,y))-1-exp(-phi.hat[3]*abs(pmax(x,y)-pmin(x,y))))
)
Vi <- phi.hat[1] * Ji + Ri + phi.hat[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){
pow <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- outer(Timei, Timei, function(x,y) exp(-abs(x - y)^pow/phi.hat[3]))  
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[4] * Ii
}

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){
kappa <- as.numeric(strsplit(process, "-")[[1]][3])
Ri <- matern(abs(outer(Timei, Timei, "-")), phi.hat[3], kappa = kappa)  
Vi <- phi.hat[1] * Ji + phi.hat[2] * Ri + phi.hat[4] * Ii
}

sum.llik <- sum.llik + dmvt(as.numeric(Yi), delta = as.numeric(Xi %*% alpha.hat), sigma = Vi, df = dof.hat, log = TRUE)

}#for(i in idlist)


### preparing the output

output <- list()

if(process == "bm"){

est.em <- cbind(param.new, c(sqrt(diag(solve(hessian))), NA, NA, NA, NA), sqrt(diag(solve(-hessian.full))))
est.em.z <- est.em[,1]/est.em[,2]
est.em.p <- pnorm(-abs(est.em.z))*2
est.em   <- cbind(est.em, est.em.z, est.em.p)
rownames(est.em) <- c(colnames(x.em), "omegasq", "sigmasq", "tausq", "dof")
colnames(est.em) <- c("Estimate", "SE-theoretical", "SE-cd", "Z", "p")

output$title      <- "Linear mixed-effects model with t-distribution and Brownian motion"
output$est        <- est.em
output$max.loglik <- sum.llik
output$varcov.alpha.theoretical <- solve(hessian)
output$varcov.full.central.difference <- solve(-hessian.full) 
rownames(output$varcov.full.central.difference) <- 
  colnames(output$varcov.full.central.difference) <- c(colnames(hessian), "omegasq", "sigmasq", "tausq", "dof")


}#bm

if(process == "ibm"){

est.em <- cbind(param.new, c(sqrt(diag(solve(hessian))), NA, NA, NA, NA), sqrt(diag(solve(-hessian.full))))
est.em.z <- est.em[,1]/est.em[,2]
est.em.p <- pnorm(-abs(est.em.z))*2
est.em   <- cbind(est.em, est.em.z, est.em.p)
rownames(est.em) <- c(colnames(x.em), "omegasq", "sigmasq", "tausq", "dof")
colnames(est.em) <- c("Estimate", "SE-theoretical", "SE-cd", "Z", "p")

output$title      <- "Linear mixed-effects model with t-distribution and integrated Brownian motion"
output$est        <- est.em
output$max.loglik <- sum.llik
output$varcov.alpha.theoretical <- solve(hessian)
output$varcov.full.central.difference <- solve(-hessian.full) 
rownames(output$varcov.full.central.difference) <- 
  colnames(output$varcov.full.central.difference) <- c(colnames(hessian), "omegasq", "sigmasq", "tausq", "dof")


}#ibm

if(process == "iou"){

est.em <- cbind(param.new, c(sqrt(diag(solve(hessian))), NA, NA, NA, NA, NA), sqrt(diag(solve(-hessian.full))))
est.em.z <- est.em[,1]/est.em[,2]
est.em.p <- pnorm(-abs(est.em.z))*2
est.em   <- cbind(est.em, est.em.z, est.em.p)
rownames(est.em) <- c(colnames(x.em), "omegasq", "kappasq", "nu", "tausq", "dof")
colnames(est.em) <- c("Estimate", "SE-theoretical", "SE-cd", "Z", "p")

output$title      <- "Linear mixed-effects model with t-distribution and integrated Ornstein-Uhlenbeck process"
output$est        <- est.em
output$max.loglik <- sum.llik
output$varcov.alpha.theoretical <- solve(hessian)
output$varcov.full.central.difference <- solve(-hessian.full) 
rownames(output$varcov.full.central.difference) <- 
  colnames(output$varcov.full.central.difference) <- c(colnames(hessian), "omegasq", "kappasq", "nu", "tausq", "dof")

}#iou

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "powered"){

est.em <- cbind(param.new, c(sqrt(diag(solve(hessian))), NA, NA, NA, NA, NA), sqrt(diag(solve(-hessian.full))))
est.em.z <- est.em[,1]/est.em[,2]
est.em.p <- pnorm(-abs(est.em.z))*2
est.em   <- cbind(est.em, est.em.z, est.em.p)
rownames(est.em) <- c(colnames(x.em), "omegasq", "sigmasq", "phi", "tausq", "dof")
colnames(est.em) <- c("Estimate", "SE-theoretical", "SE-cd", "Z", "p")

output$title      <- "Linear mixed-effects model with t-distribution and powered correlation function"
output$est        <- est.em
output$max.loglik <- sum.llik
output$varcov.alpha.theoretical <- solve(hessian)
output$varcov.full.central.difference <- solve(-hessian.full) 
rownames(output$varcov.full.central.difference) <- 
  colnames(output$varcov.full.central.difference) <- c(colnames(hessian), "omegasq", "sigmasq", "phi", "tausq", "dof")


}#powered

if(strsplit(process, "-")[[1]][1] == "sgp" & strsplit(process, "-")[[1]][2] == "matern"){

est.em <- cbind(param.new, c(sqrt(diag(solve(hessian))), NA, NA, NA, NA, NA), sqrt(diag(solve(-hessian.full))))
est.em.z <- est.em[,1]/est.em[,2]
est.em.p <- pnorm(-abs(est.em.z))*2
est.em   <- cbind(est.em, est.em.z, est.em.p)
rownames(est.em) <- c(colnames(x.em), "omegasq", "sigmasq", "phi", "tausq", "dof")
colnames(est.em) <- c("Estimate", "SE-theoretical", "SE-cd", "Z", "p")

output$title      <- "Linear mixed-effects model with t-distribution and Matern correlation function"
output$est        <- est.em
output$max.loglik <- sum.llik
output$varcov.alpha.theoretical <- solve(hessian)
output$varcov.full.central.difference <- solve(-hessian.full) 
rownames(output$varcov.full.central.difference) <- 
  colnames(output$varcov.full.central.difference) <- c(colnames(hessian), "omegasq", "sigmasq", "phi", "tausq", "dof")

}#matern

output

}


