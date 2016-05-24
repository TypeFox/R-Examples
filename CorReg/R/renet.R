# ' correction for enet not working when X has only one variable
renet<-function (x, y, lambda = 0, max.steps, normalize = TRUE, intercept = TRUE, 
          trace = FALSE, eps = .Machine$double.eps) 
{
call <- match.call()
nm <- dim(x)
n <- nm[1]
m <- nm[2]
im <- seq(m)
one <- rep(1, n)
vn <- dimnames(x)[[2]]
meanx <- drop(one %*% x)/n
if (intercept == FALSE) {
   meanx <- rep(0, m)
}
x <- scale(x, meanx, FALSE)
normx <- sqrt(drop(one %*% (x^2)))
if (normalize == FALSE) {
   normx <- rep(1, m)
}
if (any(normx < eps * sqrt(n))) 
   stop("Some of the columns of x have zero variance")
names(normx) <- NULL
x <- scale(x, FALSE, normx)
mu <- mean(y)
if (intercept == FALSE) {
   mu <- 0
}
y <- drop(y - mu)
d1 <- sqrt(lambda)
d2 <- 1/sqrt(1 + lambda)
Cvec <- drop(t(y) %*% x) * d2
ssy <- sum(y^2)
residuals <- c(y, rep(0, m))
if (lambda > 0) {
   maxvars <- m
}
if (lambda == 0) {
   maxvars <- min(m, n - 1)
}
if (missing(max.steps)) {
   max.steps <- 50 * maxvars
}
L1norm <- 0
penalty <- max(abs(Cvec))
beta <- rep(0, m)
betactive <- list(NULL)
first.in <- integer(m)
active <- NULL
Actset <- list(NULL)
df <- 0
if (lambda != 0) {
   Cp <- ssy
}
ignores <- NULL
actions <- as.list(seq(max.steps))
drops <- FALSE
Sign <- NULL
R <- NULL
k <- 0
while ((k < max.steps) & (length(active) < maxvars)) {
   action <- NULL
   k <- k + 1
   inactive <- if (k == 1) {
      im
   }else {im[-c(active, ignores)]}
   C <- Cvec[inactive]
   Cmax <- max(abs(C))
   if (!any(drops)) {
      new <- abs(C) == Cmax
      C <- C[!new]
      new <- inactive[new]
      for (inew in new) {
         R <- updateRR(x[, inew], R, x[, active], lambda)
         if (attr(R, "rank") == length(active)) {
            nR <- seq(length(active))
            R <- R[nR, nR, drop = FALSE]
            attr(R, "rank") <- length(active)
            ignores <- c(ignores, inew)
            action <- c(action, -inew)
            if (trace) 
               cat("LARS-EN Step", k, ":\t Variable", inew, 
                   "\tcollinear; dropped for good\n")
         }else {
            if (first.in[inew] == 0) 
               first.in[inew] <- k
            active <- c(active, inew)
            Sign <- c(Sign, sign(Cvec[inew]))
            action <- c(action, inew)
            if (trace) 
               cat("LARS-EN Step", k, ":\t Variable", inew, 
                   "\tadded\n")
         }
      }
   }else action <- -dropid
   Gi1 <- backsolve(R, backsolvet(R, Sign))
   A <- 1/sqrt(sum(Gi1 * Sign))
   w <- A * Gi1
   u1 <- drop(x[, active, drop = FALSE] %*% w * d2)
   u2 <- rep(0, m)
   u2[active] <- d1 * d2 * w
   u <- c(u1, u2)
   if (lambda > 0) {
      maxvars <- m - length(ignores)
   }
   if (lambda == 0) {
      maxvars <- min(m - length(ignores), n - 1)
   }
   if (length(active) >= maxvars) {
      gamhat <- Cmax/A
   } else {
      a <- (drop(u1 %*% x[, -c(active, ignores)]) + d1 * 
               u2[-c(active, ignores)]) * d2
      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
      gamhat <- min(gam[gam > eps], Cmax/A)
      Cdrop <- c(C - gamhat * a, -C + gamhat * a) - (Cmax - 
                                                        gamhat * A)
   }
   dropid <- NULL
   b1 <- beta[active]
   z1 <- -b1/w
   zmin <- min(z1[z1 > eps], gamhat)
   if (zmin < gamhat) {
      gamhat <- zmin
      drops <- z1 == zmin
   }else drops <- FALSE
   beta[active] <- beta[active] + gamhat * w
   betactive[[k]] <- beta[active]
   Actset[[k]] <- active
   residuals <- residuals - (gamhat * u)
   Cvec <- (drop(t(residuals[1:n]) %*% x) + d1 * residuals[-(1:n)]) * 
      d2
   L1norm <- c(L1norm, sum(abs(beta[active]))/d2)
   penalty <- c(penalty, penalty[k] - abs(gamhat * A))
   if (any(drops)) {
      dropid <- seq(drops)[drops]
      for (id in rev(dropid)) {
         if (trace) 
            cat("LARS-EN Step", k, ":\t Variable", active[id], 
                "\tdropped\n")
         R <- downdateR(R, id)
      }
      dropid <- active[drops]
      beta[dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
   }
   if (!is.null(vn)) 
      names(action) <- vn[abs(action)]
   actions[[k]] <- action
}
allset <- Actset[[1]]
if(k>=2){
   for (i in 2:k) {
      allset <- union(allset, Actset[[i]])
   }
}


allset <- sort(allset)
max.p <- length(allset)
beta.pure <- matrix(0, k + 1, max.p)
for (i in 2:(k + 1)) {
   for (j in 1:length(Actset[[i - 1]])) {
      l <- c(1:max.p)[allset == Actset[[i - 1]][j]]
      beta.pure[i, l] <- betactive[[i - 1]][j]
   }
}
beta.pure <- beta.pure/d2
dimnames(beta.pure) <- list(paste(0:k), vn[allset])
k <- dim(beta.pure)[1]
df <- 1:k
for (i in 1:k) {
   a <- drop(beta.pure[i, ])
   df[i] <- 1 + length(a[a != 0])
}
residuals <- y - x[, allset, drop = FALSE] %*% t(beta.pure)
beta.pure <- scale(beta.pure, FALSE, normx[allset])
RSS <- apply(residuals^2, 2, sum)
R2 <- 1 - RSS/RSS[1]
Cp <- ((n - m - 1) * RSS)/rev(RSS)[1] - n + 2 * df
object <- list(call = call, actions = actions[seq(k)], allset = allset, 
               beta.pure = beta.pure, vn = vn, mu = mu, normx = normx[allset], 
               meanx = meanx[allset], lambda = lambda, L1norm = L1norm, 
               penalty = penalty * 2/d2, df = df, Cp = Cp, sigma2 = rev(RSS)[1]/(n - 
                                                                                    m - 1))
class(object) <- "enet"
object
}