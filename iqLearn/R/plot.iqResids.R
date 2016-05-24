plot.iqResids <-
function (x, ...){
  qqnorm (x, ...);
  qqline (x);
}
