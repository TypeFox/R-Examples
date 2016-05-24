library(freqdom)

set.seed(4)

precision = 0.0001

OP = array(0,c(3,2,2))
OP[1,,] = diag(2)
OP[2,,] = diag(2)
OP[3,,] = diag(2)

A = timedom(OP,c(-1,0,2))

X = matrix(rnorm(200),100,2)
Xt = timedom(X)

Xf = fourier.transform(Xt)
Af = fourier.transform(A)

Y = A %c% X    # time domain convolution
Cf = Af %*% Xf # frequency domain product

Yt = invfourier(Cf,1:100)
Yr = Yt$operators[,,1]

par(mfrow=c(1,1))
plot(Y[,1])
lines(Yr[,1])
