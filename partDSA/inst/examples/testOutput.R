library(partDSA)

# 1) Continuous outcome with categorical and continuous predictor variables

x <- data.frame(matrix(rnorm(200), 50, 4) )
x[,2] <- as.numeric(x[,2] > 0)
x[sample(1:50, 10), 2] <- 2
y <- x[,1]^2 / 4 + x[,2]^2 + rnorm(50, 0, .25)
x[,3] <- as.factor(x[,3] > 0)
x[,2] <- as.factor(x[,2])
x[,2] <- as.factor(ifelse(x[,2]==0, "a", ifelse(x[,2]==1, "b", "c")))
x[sample(1:50, 10), 1] <- NA
x[sample(1:50, 10), 2] <- NA
x[sample(1:50, 10), 3] <- NA
x[sample(1:50, 10), 4] <- NA
model <- partDSA(data.frame(x), y, control=DSA.control(missing="impute.at.split"))
print(model)
showDSA(model)

# 2) Categorical outcome with continuous predictor variable

x <- matrix(rnorm(300), 100, 3)
y <- as.factor(sample(c("a", "b", "c"), 100, TRUE) )
x.test <- matrix(rnorm(150), 50, 3)
y.test <- as.factor(sample(c("a", "b", "c"), 50, TRUE))
model <- partDSA(x, y)
print(model)
showDSA(model)

# 3) Continuous outcome with continuous predictor variables

x <- data.frame(matrix(rnorm(150), 50, 3))
y <- x[,1]^2 + x[,2]^2 + rnorm(50)
names(x) <- c("A", "B", "C")
model <- partDSA(x, y)
print(model)
showDSA(model)
