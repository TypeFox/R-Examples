stopifnot(require(glorius))

data(dat30)

M <- solarMultivar(traits = c("trait1", "trait2"), data = dat30)

gdat <- genocovdat30
stopifnot(all(M$ids == rownames(gdat)))

gdat <- rbind(gdat, gdat)
rownames(gdat) <- M$data$tid

A1 <- matrixAssoc(as.formula(M$formula), M$data, gdat, correlation = as.matrix(M$correlation), 
  pvOutputThreshold = 1, id.var = "tid")
head(A1$snpf)  

A2 <- matrixAssoc(as.formula(M$formula), M$data, gdat, correlation = as.matrix(M$correlation), 
  pvOutputThreshold = 1, id.var = "tid", interaction = "tnum")
head(A2$snpf)

