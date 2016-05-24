emma <- function(in.name, nlev, lower, upper, out.name, opt = "mn", nd = 10, na = 5, 
weight, C = 20, w1 = 0.7, w2 = 0.4, c1i = 2.5, c1f = 0.5, c2i = 0.5, c2f = 2.5, b = 5, 
pr.mut, graph, fn1 = NULL, fn2 = NULL, fn3 = NULL, fn4 = NULL, nresp = 1) 
{

  ## identification of initial set of experimental runs (initialization)
  tn <- emmat0(in.name, nlev, lower, upper, out.name, nd, fn1, fn2, fn3, fn4)

  ## iterative identification of additional sets of experimental runs
  for(t in 1:(C-1))  {
    tn <- emmatn(t, tn, na, opt, weight, C, w1, w2, c1i, c1f, c2i, c2f, b, pr.mut, graph, fn1, fn2, fn3, fn4, nresp)
    tn <- emmacheck(tn, graph, fn1, fn2, fn3, fn4, nresp)
  }

  return(tn)
}
