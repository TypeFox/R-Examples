## library(appell)

## compare the results of both algorithms
## for random test data.

## todo: add better tests trying to replicate published results? 

nTest <- 100L
set.seed(38)

a <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))
b <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))
c <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))
z <- complex(real=rnorm(nTest),
             imaginary=rnorm(nTest))

tableHyp2f1 <- matrix(nrow=nTest,
                      ncol=2L,
                      dimnames=
                      list(NULL,
                           c("forrey", "michel.stoitsov")))

for(i in seq_len(nTest))
{
    tableHyp2f1[i, "forrey"] <- hyp2f1(a[i], b[i], c[i], z[i],
                                       algorithm="forrey")
    tableHyp2f1[i, "michel.stoitsov"] <- hyp2f1(a[i], b[i], c[i], z[i],
                                                algorithm="michel.stoitsov")
}

tableHyp2f1

abs(tableHyp2f1[, "forrey"] - tableHyp2f1[, "michel.stoitsov"])
## so very small differences,
## at least in this range of function parameters.
