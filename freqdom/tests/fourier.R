library(freqdom)
OP = array(0,c(3,1,1))
OP[1,,] = 1
OP[3,,] = 1
A = timedom(OP,-1:1)
FT = fourier.transform(A) #should give cosine
Ai = invfourier(FT,lags=-1:1)

if (abs(sum((Ai$operators[1:3,,] - A$operators[1:3,,])^2)) > 0.001)
  stop("ERROR")