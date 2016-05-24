`propdiff.freq0` <-
function(len, c1, d1, c2, d2, level=0.95)
{
  propdiff.freq(len, c1/(c1+d1), c2/(c2+d2), level)
}

