
library(rbenchmark)

N <- 1e5
v <- vector(mode="numeric", length=N)

## res <- benchmark(zrnormMT(N),       # Marsgalia and Tsang, JSS, 2000
##                  zrnormLZLLV(N),    # Leong, Zhang et al, JSS, 2005
##                  zrnormV1(N),       # based on initial Burkardt implementation
##                  zrnormVecV1(v),    # fill a pre-supplied vector
##                  #zrnormStlV1(N),   # fill STL vector
##                  rnorm(N),          # R as a baseline
##                  zrnorm(N),         # based on updated Burkardt implementation
##                  zrnormVec(v),      # fill a pre-supplied vector
##                  #zrnormStl(N),     # fill STL vector
##                  zrnormgsl(N),      # GSL's ziggurat by Voss
##                  zrnormV1b(N),      # based on initial Burkardt impl, mod'ed
##                  replications=1000, order="relative")
## print(res[,1:4])

res <- benchmark(zrnormMT(N),       # Marsgalia and Tsang, JSS, 2000
                 zrnormLZLLV(N),    # Leong, Zhang et al, JSS, 2005
                 #zrnormV1(N),       # based on initial Burkardt implementation
                 #zrnormVecV1(v),    # fill a pre-supplied vector
                 #zrnormStlV1(N),   # fill STL vector
                 #rnorm(N),          # R as a baseline
                 zrnorm(N),         # based on updated Burkardt implementation
                 zrnormVec(v),      # fill a pre-supplied vector
                 #zrnormStl(N),     # fill STL vector
                 zrnormGSL(N),      # GSL's ziggurat by Voss
                 #zrnormV1b(N),      # based on initial Burkardt impl, mod'ed
                 zrnormQL(N),       # QuantLib variant
                 zrnormGl(N),       # Gretl
                 replications=1000, order="relative")
print(res[,1:4])

if (requireNamespace("microbenchmark", quietly=TRUE)) {
    res <- microbenchmark(zrnormMT(N), zrnorm(N), zrnormLZLLV(N), zrnormGSL(N), zrnormQL(N),
                          zrnormGl(N), zrnormV1(N), zrnormV1b(N), rnorm(N),
                          times=1000, control=list(warmup=20))
    oo <- order(summary(res)[,"median"])
    res$expr <- ordered(x=as.numeric(res$expr),
                        levels=oo,
                        labels=levels(res$expr)[oo])
    print(res)
    if (interactive())
        if (requireNamespace("ggplot2", quietly=TRUE))
            ggplot2::autoplot(res)
}
