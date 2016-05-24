phs <- cbind(c(104, 189), c(10933, 10845)); phs
rownames(phs) <- c("aspirin", "placebo")                  # add row names
colnames(phs) <- c("heart attack", "no heart attack")     # add column names
phs 
xchisq.test(phs)

