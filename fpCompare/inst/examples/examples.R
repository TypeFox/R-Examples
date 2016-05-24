x1 <- 0.5 - 0.3
x2 <- 0.3 - 0.1
x1 == x2                           # FALSE on most machines
x1 %==% x2                         # TRUE everywhere
identical(all.equal(x1, x2), TRUE) # TRUE everywhere

set.seed(123)
a <- 1:6
b <- jitter(1:6, 1e-7)
print(rbind(a,b), digits=16)

b %<=% a
b %<<% a
b %>=% a
b %>>% a
b %==% a
b %!=% a
