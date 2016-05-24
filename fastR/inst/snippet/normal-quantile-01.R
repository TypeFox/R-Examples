x <- c(-0.16,1.17,-0.43,-0.02,1.06,-1.35,0.65,-1.12,0.03,-1.44)
# sort the data
x.sorted <- sort(x); x.sorted
q <- seq(0.05,0.95, by=0.1); q
y <- qnorm(q); y
