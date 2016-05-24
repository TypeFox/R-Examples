`mu.freq` <-
function(len, lambda, level=0.95)
{
  ceiling(4/lambda/len/len*qnorm((1+level)/2)^2)
}
