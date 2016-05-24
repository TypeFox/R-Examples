boot.nm <- function(formula, data, id, timeVar, result, matern = TRUE, kappa.or.power, nboot = 100, 
                        tol.lmenssp = 1e-8, maxiter.lmenssp = 500){

mf <- model.frame(formula = formula, data = data)
#y  <- as.matrix(model.extract(mf, "response"))
x   <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))

ncov <- ncol(x)

alpha   <- result[1 : ncov,   1, drop = F]
omegasq <- result[(ncov + 1), 1, drop = T]
sigmasq <- result[(ncov + 2), 1, drop = T]
phi     <- result[(ncov + 3), 1, drop = T]
tausq   <- result[(ncov + 4), 1, drop = T]

numb.obs  <- tapply(id, id, function(x) length(x))
numb.subj <- length(unique(id)) 
id.unique <- unique(id)

save <- c()

for (boot in 1 : nboot){

cat("Bootstrap replication = ", boot, "\n")

U <- rep(rnorm(numb.subj, 0, sqrt(omegasq)), numb.obs)

W <- c()
 
time.id <- data.frame(timeVar = timeVar, id = id)

for (i in 1 : numb.subj){

fui <- time.id[time.id[, "id"] == id.unique[i], "timeVar"]
if(matern == TRUE) Ri  <- sigmasq * matern(abs(outer(fui, fui, "-")), phi = phi, kappa = kappa.or.power)
if(matern == FALSE) Ri <- sigmasq * outer(fui, fui, function(x,y) exp(-abs(x - y)^kappa.or.power/phi)) 
Wi  <- rmvnorm(1, rep(0, ncol(Ri)), Ri)
W   <- c(W, Wi)

} 

Z <- rnorm(length(W), 0, sqrt(tausq))

data$resp     <- x %*% alpha + U + W + Z

if(length(strsplit(as.character(formula), "~")) == 2){
formula2 <- as.formula(paste("resp", strsplit(as.character(formula), "~")[[2]], sep = "~"))
}
if(length(strsplit(as.character(formula), "~")) == 3){
formula2 <- as.formula(paste("resp", strsplit(as.character(formula), "~")[[3]], sep = "~"))
}

if(matern == TRUE) process  <- paste("sgp-matern", kappa.or.power, sep = "-")
if(matern == FALSE) process <- paste(paste("sgp-matern", kappa.or.power, sep = "-"), "nm", sep = "-") 

init     <- c(log(omegasq/sigmasq), log(phi), log(tausq/sigmasq))

fit.matern <- lmenssp(formula = formula2, data = data, id = id, timeVar = timeVar, 
                      process = process, init = init, tol = tol.lmenssp, maxiter = maxiter.lmenssp)

save <- rbind(save, fit.matern$estimate[, 1])

}#for (boot in 1 : nboot)

new.result <- cbind(result[, 1, drop = F], colMeans(save), result[, 2, drop = F], 
                    sqrt(diag(cov(save))), result[, 3 : 4, drop = F])
colnames(new.result)[2 : 6] <- c("Bootstrap Mean", "Theoretical SE", 
                                 "Bootstrap SE", "Theoretical Z", "Theoretical p")
 
output <- list()
if(matern == TRUE) output$title <- "Boostrapping for mixed model with stationary process - Matern correlation"
if(matern == FALSE) output$title <- "Boostrapping for mixed model with stationary process - powered correlation"
output$result <- new.result
output$varcov <- cov(save)
output

}#boot.sgp

