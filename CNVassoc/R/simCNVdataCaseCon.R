simCNVdataCaseCon <- function (n0,n1, w0, or, mu.surrog0, sd.surrog0, mu.surrog1=mu.surrog0, sd.surrog1=sd.surrog0, random=TRUE)
{
    k <- length(w0)
    w0 <- w0/sum(w0)
    or <- c(1,or)
    w1 = or*w0/sum(or*w0)
    r0 <- rmultinom(1,n0,w0)
    r1 <- rmultinom(1,n1,w1)
    s0 <- unlist(sapply(seq(along=r0),function(j) rnorm(r0[j], mean=mu.surrog0[j], sd = sd.surrog0[j])))
    s1 <- unlist(sapply(seq(along=r1),function(j) rnorm(r1[j], mean=mu.surrog1[j], sd = sd.surrog1[j])))
    cnv0 <- rep(1:k, r0)
    cnv1 <- rep(1:k, r1)
    out <- data.frame(resp=rep(0:1,c(n0,n1)), cnv=c(cnv0,cnv1), surrog = c(s0,s1))
    if (random) out <- out[sample(1:nrow(out)), ]
    out
}