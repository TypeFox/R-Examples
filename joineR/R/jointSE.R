jointSE <- function(fitted, n.boot, gpt, lgpt, max.it, tol, print.detail = FALSE){

data <- fitted$data
id <- fitted$data$subj.col
time.long <- fitted$data$time.col
q <- length(diag(fitted$sigma.u))
paranames <- c(row.names(fitted$coefficients$fixed$longitudinal), names(fitted$coefficients$fixed$survival),
    names(fitted$coefficients$latent), paste("U_", 0:(q-1), sep = "") , "Residual")
compnames <- rep("", length(paranames))
compnames[1] <- "Longitudinal"
lb1 <- length(fitted$coefficients$fixed$longitudinal[,1])
lb2 <- length(fitted$coefficients$fixed$survival)
lg <- length(fitted$coefficients$latent)
compnames[lb1 + 1] <- "Survival"
compnames[lb1 + lb2 + 1] <- "Association"
compnames[lb1 + lb2 + lg + 1] <- "Variance"
if (missing(gpt)) {
    gpt <- 3
}
if (missing(lgpt)) {
    lgpt <- 10
}
if (missing(max.it)) {
    max.it <- 200
}
if (missing(tol)) {
    tol <- 0.001
}
model <- fitted$model
surv.formula <- fitted$formulae$sformula
long.formula <- fitted$formulae$lformula
sepassoc <- fitted$sepassoc
data.surv <- cbind(fitted$data$survival, fitted$data$baseline)
surv.frame <- model.frame(surv.formula, data = data.surv)
if (dim(surv.frame)[2] == 1){
n.est <- dim(as.matrix(fitted$coefficients$fixed$longitudinal))[1] +
dim(as.matrix(fitted$coefficients$latent))[1] +
dim(as.matrix(diag(fitted$sigma.u)))[1] + 1
} else  { 
n.est <- dim(as.matrix(fitted$coefficients$fixed$longitudinal))[1] +
dim(as.matrix(fitted$coefficients$fixed$survival))[1] +
dim(as.matrix(fitted$coefficients$latent))[1] +
dim(as.matrix(diag(fitted$sigma.u)))[1] + 1
}
out <- matrix(0, n.boot + 2, n.est)
nsubj <- length(fitted$data$subject)
for (i in 1:n.boot) {
    s.new <- sample.jointdata(data, nsubj, replace = TRUE)
    fitb <- joint(data = s.new, long.formula = long.formula, surv.formula = surv.formula, model = model, sepassoc = sepassoc, 
              gpt = gpt, max.it = max.it, tol = tol, lgpt = lgpt)
    b1 <- as.numeric(as.vector(as.matrix(fitb$coefficients$fixed$longitudinal[,1])))
    b3 <- as.numeric(as.vector(as.matrix(fitb$coefficients$latent)))
    b4 <- as.numeric(as.vector(as.matrix(diag(fitb$sigma.u))))
    b5 <- as.numeric(as.vector(as.matrix(fitb$sigma.z)))
    if (dim(surv.frame)[2] != 1) { 
        b2 <- as.numeric(as.vector(as.matrix(fitb$coefficients$fixed$survival)))
        out[i, ] <- c(b1, b2, b3, b4, b5)
        ests <- out[i, ]
        if (print.detail) { 
            detail <- data.frame(iteration = i, t(ests))
            names(detail) <- c("Iteration", paranames)
            print(detail)
        }
    } 
    else {
        out[i, ] <- c(b1, b3, b4, b5)
        ests <- out[i, ]
        if (print.detail) {
            detail <- data.frame(iteration = i, t(ests))
            names(detail) <- c("Iteration", paranames)
	    print(detail)
        }
    }
}
    i <- 1
    while (out[i, 1] != 0) i = i + 1
    out <- out[1 : (i - 1), ]
    se <- 0
    ci1 <- 0
    ci2 <- 0
    if (n.boot == 1) {out <- matrix(out, nrow = 1)}
    for (i in 1:length(out[1, ])) {
        se[i] <- sqrt(var(as.numeric(out[, i])))
        if (n.boot < 100) {
       	    ci1[i] <- 0
            ci2[i] <- 0
        }
        else {
            ci1[i] <- sort(as.numeric(out[, i]))[2.5/100 * n.boot]
            ci2[i] <- sort(as.numeric(out[, i]))[97.5/100 * n.boot]
        }
    }
if (dim(surv.frame)[2] != 1){
b1 <- data.frame(cbind(compnames, paranames,
round(c(as.numeric(as.vector(as.matrix(fitted$coefficients$fixed$longitudinal))),
as.numeric(as.vector(as.matrix(fitted$coefficients$fixed$survival))),
as.numeric(as.vector(as.matrix(fitted$coefficients$latent))),
as.numeric(as.vector(as.matrix(diag(fitted$sigma.u)))),
as.numeric(as.vector(as.matrix(fitted$sigma.z)))), 4),
round(cbind(se), 4), round(ci1, 4), round(ci2, 4)))
} else {
b1 <- data.frame(cbind(compnames, paranames,
round(c(as.numeric(as.vector(as.matrix(fitted$coefficients$fixed$longitudinal))),
as.numeric(as.vector(as.matrix(fitted$coefficients$latent))),
as.numeric(as.vector(as.matrix(fitted$sigma.z))),
as.numeric(as.vector(as.matrix(diag(fitted$sigma.u))))), 4),
round(cbind(se), 4), round(ci1, 4), round(ci2, 4)))
}
    names(b1)[1:6] <- c("Component", "Parameter", "Estimate", "SE", "95%Lower", "95%Upper")
    b1
}
