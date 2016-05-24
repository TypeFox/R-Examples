# A. Looping on objects other than vectors
set.seed(12); da <- data.frame(y = rnorm(50), x = rnorm(50))
set.seed(34); db <- data.frame(y = rnorm(50), x = rnorm(50))

# A1. "list" approach
data.all <- list(da = da, db = db)
reg1 <- list() 
coe1 <- data.frame(matrix(data = 0, nrow = length(data.all), ncol = 2))
for (j in 1:length(data.all)) {
  fit <- lm(formula = y ~ x, data = data.all[[j]] )
  reg1[[j]] <- fit
  coe1[j, ] <- coef(fit)
}
coe1

# A2. "get" approach
identical(get("da"), da)  # TRUE
data.name <- c("da", "db")      
reg2 <- list()
coe2 <- NULL
for (k in data.name) {  
  w <- get(k)
  fit <- lm(formula = y ~ x, data = w)
  reg2[[k]] <- fit
  coe2 <- rbind(coe2, coef(fit))
}
coe2  

# B. Preallocation versus redimension
input <- matrix(data = 1:12, nrow = 3, ncol = 4)
nam <- c("outA", "outB", "outC", "outD")
out.prea <- data.frame(matrix(data = 0, nrow = 3, ncol = length(nam)))
colnames(out.prea) <- nam
out.redi <- temp <- NULL

for (n in 1:4) {  
  for (m in 1:3) {
    out.prea[m, n] <- input[m, n] * m * n
    temp <- c(temp, input[m, n] * m * n)  # a growing vector
  } 
  out.redi <- data.frame(matrix(data = temp, nrow = 3)) # a growing matrix
} 
out.prea
out.redi