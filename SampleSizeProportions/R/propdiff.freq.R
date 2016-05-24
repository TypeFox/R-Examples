`propdiff.freq` <-
function(len, p1.estimate, p2.estimate, level=0.95)
{
  ceiling(4*qnorm((1+level)/2)^2*(p1.estimate*(1-p1.estimate) + p2.estimate*(1-p2.estimate))/len^2)
}

