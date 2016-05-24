`ddst.uniform.test` <-
function(x, base = ddst.base.legendre, c = 2.4, B=1000, compute.p = F, Dmax = 10, ...) {

# method.name = as.character(substitute(base)) 
 # only Legendre is implemented yet
 base = ddst.base.legendre
 method.name = "ddst.base.legendre"

  n = length(x)
if(n<5) 
	stop("length(x) should be at least 5")
  coord = ddst.uniform.Nk(x, base, Dmax = Dmax)    # coord square times n
  l = ddst.IIC(coord, n, c)
  attr(l, "names") = "n. coord"
  t = coord[l]
  attr(t, "names") = "WT"
  result = list(statistic = t, parameter = l, method = "Data Driven Smooth Test for Uniformity")
  result$data.name = paste(paste(as.character(substitute(x)), collapse=""), ",   base: ", method.name, "   c: ", c, sep="")
  class(result) = "htest"
  if (compute.p) {
     tmp = numeric(B)
     for (i in 1:B) {
        y = runif(n)
        tmpC = ddst.uniform.Nk(y, base, Dmax = Dmax)
        l = ddst.IIC(tmpC, n, c)
        tmp[i] = tmpC[l]
     }
     result$p.value = mean(tmp > t)
  }
  result
}

