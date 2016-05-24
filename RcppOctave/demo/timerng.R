
library(RcppOctave)
require(rbenchmark)			# to benchmark examples

## define a Gamma RNG draw function
o_source(text="
  function x = orndg(a, b, n)
    x = b * randg(a, n, 1);
  end
")

## define a Normal RNG draw function
o_source(text="
  function x = orndn(n)
    x = randn(n,1);
  end
")

N <- 500

set.seed(42)  # reset RNG
x1 <- c(.O$orndg(9,1,N))
set.seed(42)  # reset RNG
y1 <- rgamma(N,9,1)
stopifnot(all.equal(x1, y1))

res <- benchmark(.O$orndg(9,1,N),
                 rgamma(N,9,1),
                 o_rgamma(9,N,1),
                 columns = c("test", "replications", "elapsed", "relative"),
                 order="relative",
                 replications=1000)
print(res)

set.seed(42)  # reset RNG
x1 <- c(.O$orndn(N,1))
set.seed(42)  # reset RNG
y1 <- rnorm(N)
stopifnot(all.equal(x1, y1))

res <- benchmark(.O$orndn(N),
                 rnorm(N),
                 o_rnorm(N,1),
                 columns = c("test", "replications", "elapsed", "relative"),
                 order="relative",
                 replications=1000)
print(res)
