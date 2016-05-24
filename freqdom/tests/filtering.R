library(freqdom)

set.seed(4)

precision = 0.0001

# Tes 1t: X_t = Y_{t-1} + Y_{t} + Y_{t+1}
OP = array(0,c(3,2,2))
OP[1,,] = diag(2)
OP[2,,] = diag(2)
OP[3,,] = diag(2)
A = timedom(OP,-1:1)

X = matrix(rnorm(200),100,2)
Y = filter.process(X,A)

if (abs(Y[4,2] - X[3,2] - X[4,2] - X[5,2]) > precision)
  stop("error test 1")

# Test: X_t = 2*Y_{t-1} + Y_{t} + Y_{t+1}
A$operators[3,,] = diag(2)*2
A$lags = c(-1,0,1)

Y = filter.process(X,A)

if (abs(Y[4,2] - 2*X[3,2] - X[4,2] - X[5,2]) > precision)
  stop("error")

# Test: X_t = 2*Y_{t-2} + Y_{t} + Y_{t+1}
A$operators[3,,] = diag(2)*2
A$lags = c(-1,0,2)

Y = filter.process(X,A)

if (abs(Y[4,2] - 2*X[2,2] - X[4,2] - X[5,2]) > precision)
  stop("error")

# Test: X_ti = 2*Y_{t-2}{2-i} + Y_ti + Y_{t+1}i
A$operators[3,,] = 2*matrix(c(0,1,1,0),2,2)
A$lags = c(-1,0,2)

Y = filter.process(X,A)

if (abs(Y[4,2] - 2*X[2,1] - X[4,2] - X[5,2]) > precision)
  stop("error")

# Test: X_ti = 2*Y_{t-2000}{2-i} + Y_ti + Y_{t+1}i
A$lags = c(-1,0,2000)

Y = filter.process(X,A)

if (abs(Y[4,2] - X[4,2] - X[5,2]) > precision)
  stop("error")
