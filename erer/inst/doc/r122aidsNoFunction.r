# 0. load library; inputs and choices
library(systemfit); library(erer) 
data(daBedRaw); colnames(daBedRaw)
wa <- c(2001, 1); wb <- c(2008, 12); choice <- c("CN", "VN", "ID")

# 1. With two user-defined functions: aiData(); aiStaFit()
pit <- aiData(x = daBedRaw, label = choice, start = wa, end = wb)
cow <- pit$out; round(head(cow), 3)

sh <- paste("s",   c(choice, "RW"), sep = "")
pr <- paste("lnp", c(choice, "RW"), sep = "")
rr <- aiStaFit(y = cow, share = sh, price = pr, expen = "rte", 
  hom = TRUE, sym = TRUE)
summary(rr)
names(rr); rr$formula
rr$res.matrix

# 2. Without two user-defined functions: aiData(); aiStaFit()
# 2.1 Prepare data for AIDS
vn2 <- paste("v", choice, sep = "")
qn2 <- paste("q", choice, sep = "")
x <- window(daBedRaw, start = wa, end = wb)
y <- x[, c(vn2, "vWD", qn2, "qWD")]
vRW <- y[, "vWD"] - rowSums(y[, vn2])
qRW <- y[, "qWD"] - rowSums(y[, qn2])    
value <- ts.union(y[, vn2], vRW); colnames(value) <- c(vn2, "vRW")
quant <- ts.union(y[, qn2], qRW); colnames(quant) <- c(qn2, "qRW")

price <- value / quant; colnames(price) <- c("pCN", "pVN", "pID", "pRW")
lnp <- log(price); colnames(lnp) <- c("lnpCN", "lnpVN", "lnpID", "lnpRW")
m <- ts(rowSums(value), start = wa, end = wb, frequency = 12)
share <- value / m; colnames(share) <- c("sCN", "sVN", "sID", "sRW")

rte <- log(m) - rowSums(share * lnp)
dee <- ts.union(share, rte, lnp)
colnames(dee) <- c(colnames(share), "rte", colnames(lnp))
round(head(dee), 3)
identical(cow, dee)  # TRUE

# 2.2  Formula and restriction matrix
mod <- list(China = sCN ~ 1 + rte + lnpCN + lnpVN + lnpID + lnpRW,
          Vietnam = sVN ~ 1 + rte + lnpCN + lnpVN + lnpID + lnpRW,
        Indonesia = sID ~ 1 + rte + lnpCN + lnpVN + lnpID + lnpRW)
res.left <- matrix(data = 0, nrow = 6, ncol = 18)
res.left[1, 3:6] <- res.left[2, 9:12] <- res.left[3, 15:18] <- 1
res.left[4,   4] <- res.left[5,    5] <- res.left[6,    11] <- 1
res.left[4,   9] <- res.left[5,   15] <- res.left[6,    16] <- -1
res.right <- rep(0, times = 6)
identical(res.left, rr$res.matrix)  # TRUE 
 
# 2.3 Fit AIDS model
dd <- systemfit(formula = mod, method = "SUR", data = dee, 
  restrict.matrix = res.left, restrict.rhs = res.right)
round(summary(dd, equations = FALSE)$coefficients, digits = 3)