`predict.mar1s` <-
function(object, n.ahead = 1, start.time = 0, xreg.absdata = NULL,
         init.absdata = NULL, probs = c(0.05, 0.5, 0.95),
         n.sim = 1000, ...)
{
  x <- sim.mar1s(object, n.ahead, n.sim, start.time, xreg.absdata,
                 init.absdata)
  return(quantile(x, probs = probs))
}
