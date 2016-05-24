# 1. Data preparation
library(systemfit); library(erer); data(daBedRaw)
wa <- c(2001, 1); wb <- c(2008, 12); exporter <- c("CN", "VN", "ID")

# 1.1. Data for AIDS model generated from aiData() 
pit <- aiData(x = daBedRaw, label = exporter, start = wa, end = wb)
cow <- pit$out; round(head(cow), 3)

# 1.2.  Formula and restriction matrix
mod <- list(China = sCN ~ 1 + rte + lnpCN + lnpVN + lnpID + lnpRW,
          Vietnam = sVN ~ 1 + rte + lnpCN + lnpVN + lnpID + lnpRW,
        Indonesia = sID ~ 1 + rte + lnpCN + lnpVN + lnpID + lnpRW)
res.left <- matrix(data = 0, nrow = 6, ncol = 18)
res.left[1, 3:6] <- res.left[2, 9:12] <- res.left[3, 15:18] <- 1
res.left[4, 4] <- res.left[5,  5] <- res.left[6, 11] <- 1 
res.left[4, 9] <- res.left[5, 15] <- res.left[6, 16] <- -1
res.right <- rep(0, times = 6)
 
# 2. Estimate AIDS model with systemfit()
da <- systemfit(formula = mod, method = "OLS", data = cow,
  control = systemfit.control(singleEqSigma = TRUE))
round(summary(da, equations = FALSE)$coefficients, digits = 4)  # OLS

db <- systemfit(formula = mod, method = "SUR", data = cow, 
  restrict.matrix = res.left, restrict.rhs = res.right,
  control = systemfit.control(residCovRestricted = FALSE))
round(summary(db, equations = FALSE)$coefficients, digits = 4)  # SUR

# 3. Estimate AIDS model with matrix manipulation and GLS formulas
# 3.1. stacked y and x matrices
name.y <- paste("s", exporter, sep = "")
name.x <- c("rte", paste("lnp", c(exporter, "RW"), sep = ""))
y <- matrix(data = c(cow[, name.y]), ncol = 1)
x <- diag(nrow = length(exporter)) %x% cbind(1, cow[, name.x]) # Kronecker

# 3.2. OLS regression coefficients
ols.coe <- solve(a = crossprod(x), b = crossprod(x, y))
round(c(ols.coe), digits = 4)  # My OLS 
ols.re <- y - x %*% ols.coe
re <- matrix(data = ols.re, ncol = length(exporter))

# 3.3. SUR sigma matrix
sig <- matrix(data = 0, nrow = length(exporter), ncol = length(exporter))
for (i in 1:length(exporter)) {
  for (j in 1:length(exporter)) {
    sig[i, j] <- crossprod(re[, i], re[, j]) / (nrow(cow) - 
      (length(name.x) + 1))
  }
}
sur.sig <- sig %x% diag(nrow = nrow(cow))  # Kronecker product

# 3.4. SUR coefficients
xsx <- t(x) %*% solve(sur.sig) %*% x
r0 <- matrix(data = 0, nrow = nrow(res.left), ncol = nrow(res.left))
left <- rbind(cbind(xsx, t(res.left)), cbind(res.left, r0))
righ <- matrix(data = c(t(x) %*% solve(sur.sig) %*% y, res.right), 
  ncol = 1)
sur.coe <- solve(a = left, b = righ)[1:18, ]
round(sur.coe, digits = 4)  # My SUR