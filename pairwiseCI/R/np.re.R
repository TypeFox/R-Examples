np.re <- function (x, y, conf.level = 0.95,
 alternative = c("two.sided", "less", "greater"),
 method = c("logit", "probit", "normal", "t.app")) 
{
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    conflevel <- 1-conf.level
    if (conflevel >= 1 | conflevel <= 0) { stop("The confidence level must be between 0 and 1!") }
    alternative <- match.arg(alternative)
    method <- match.arg(method)

    ssq <- function(x) {sum(x * x)}
    logit <- function(p) {log(p/(1 - p))}
    probit <- function(p) {qnorm(p)}
    expit <- function(G) {exp(G)/(1 + exp(G))}

    BF <- function(X, n1, n2) {
        N <- n1 + n2
        daten <- X
        rdaten <- rank(daten)
        rdaten1 <- rdaten[1:n1]
        rdaten2 <- rdaten[(n1 + 1):N]
        pd <- 1/n1 * (mean(rdaten2) - (n2 + 1)/2)
        pd1 <- (pd == 1)
        pd0 <- (pd == 0)
        pd[pd1] <- 0.999
        pd[pd0] <- 0.001
        rges <- rank(c(daten[1:n1]))
        rk <- rank(c(daten[(n1 + 1):N]))
        z1 <- 1/n2 * (rdaten1 - rges)
        z2 <- 1/n1 * (rdaten2 - rk)
        sqij <- var(z1)
        sqji <- var(z2)
        vd.bf <- N * sqij/n1 + N * sqji/n2
        singular.bf <- (vd.bf == 0)
        vd.bf[singular.bf] <- N/(2 * n1 * n2)
        t.bf <- sqrt(N) * (pd - 1/2)/sqrt(vd.bf)
        df.sw <- (sqij/n1 + sqji/n2)^2/((sqij/n1)^2/(n1 - 1) + 
            (sqji/n2)^2/(n2 - 1))
        df.sw[is.nan(df.sw)] <- 1000
        erg <- c(pd, vd.bf, t.bf, df.sw)
        erg
    }

fl <- c("x","y")
cmpid <- paste("p(", fl[1], ",", fl[2], ")", sep = "")

nx <- length(x)
ny <- length(y)
n <- c(nx, ny)
ntotal <- sum(n)
XY <- c(x,y)
    test <- BF(X=XY, n1=n[1], n2=n[2])
    pd <- test[1]
    vd.bf <- test[2]
    estimate <- pd

    switch(method,

 "logit" = {
    logit.pd <- logit(pd)
    logit.dev <- 1/(pd * (1 - pd))
    vd.logit <- logit.dev^2 * vd.bf
    t.logit <- (logit.pd) * sqrt(ntotal/vd.logit)

    switch(alternative,
     "two.sided"={z.bfn <- qnorm(1 - conflevel/2)
        lower.logit <- expit(logit.pd - sqrt(vd.logit/ntotal) * z.bfn)
        upper.logit <- expit(logit.pd + sqrt(vd.logit/ntotal) * z.bfn)
        p.bflogiti <- pnorm(t.logit)
        p.bflogit <- min(2 - 2 * p.bflogiti, 2 * p.bflogiti)},
     "greater"={z.bfn <- qnorm(1 - conflevel)
        lower.logit <- expit(logit.pd - sqrt(vd.logit/ntotal) * z.bfn)
        upper.logit <- 1
        p.bflogit <- 1 - pnorm(t.logit)},
     "less"={z.bfn <- qnorm(1 - conflevel)
        upper.logit <- expit(logit.pd + sqrt(vd.logit/ntotal) * z.bfn)
	lower.logit <- 0
        p.bflogit <- pnorm(t.logit)})
  parameter <- NULL
  conf.int <- c(lower.logit, upper.logit)
  t.value <- t.logit
  names(t.value)<-"logit(t)"
  p.value <- p.bflogit
  method <- "Relative Effects, Delta-Method (Logit)"
    },

 "probit" = {

    probit.pd <- qnorm(pd)
    probit.dev <- sqrt(2 * pi)/(exp(-0.5 * qnorm(pd) * qnorm(pd)))
    vd.probit <- probit.dev^2 * vd.bf
    t.probit <- (probit.pd) * sqrt(ntotal/vd.probit)

    switch(alternative,
     "two.sided"={z.bfn <- qnorm(1 - conflevel/2)
        lower.probit <- pnorm(probit.pd - sqrt(vd.probit/ntotal) * z.bfn)
        upper.probit <- pnorm(probit.pd + sqrt(vd.probit/ntotal) * z.bfn)
        p.bfprobiti <- pnorm(t.probit)
        p.bfprobit <- min(2 - 2 * p.bfprobiti, 2 * p.bfprobiti)},
    "greater"={z.bfn <- qnorm(1 - conflevel)
        lower.probit <- pnorm(probit.pd - sqrt(vd.probit/ntotal) * z.bfn)
	upper.probit <- 1
        p.bfprobit <- 1 - pnorm(t.probit)},
     "less"={ z.bfn <- qnorm(1 - conflevel)
        upper.probit <- pnorm(probit.pd + sqrt(vd.probit/ntotal) * z.bfn)
        lower.probit <- 0
        p.bfprobit <- pnorm(t.probit)})
  parameter <- NULL
  conf.int <- c(lower.probit, upper.probit)
  t.value <- t.probit
  names(t.value)<-"probit(t)"
  p.value <- p.bfprobit
  method <- "Relative Effects, Delta-Method (Probit)"
    },

 "normal" = {
    t.bf <- test[3]
 switch(alternative, 
   "two.sided"={z.bfn <- qnorm(1 - conflevel/2)
        lower.bfn <- pd - sqrt(vd.bf/ntotal) * z.bfn
        upper.bfn <- pd + sqrt(vd.bf/ntotal) * z.bfn
        p.bfni <- pnorm(t.bf)
        p.bfn <- min(2 - 2 * p.bfni, 2 * p.bfni)},
   "greater"={z.bfn <- qnorm(1 - conflevel)
        lower.bfn <- pd - sqrt(vd.bf/ntotal) * z.bfn
        upper.bfn = 1
        p.bfn <- 1 - pnorm(t.bf)},
   "less"={z.bfn <- qnorm(1 - conflevel)
        upper.bfn <- pd + sqrt(vd.bf/ntotal) * z.bfn
        lower.bfn <- 0
        p.bfn <- pnorm(t.bf)})
parameter <- NULL
conf.int <- c(lower.bfn, upper.bfn)
t.value <- t.bf
names(t.value)<-"t.bf"
p.value <- p.bfn
method <- "Relative Effects, Normal Approximation"
    },

 "t.app" = {
    t.bf <- test[3]
    df.sw <- test[4]
    switch( alternative,
     "two.sided"={z.bft <- qt(1 - conflevel/2, df.sw)
        lower.bft <- pd - sqrt(vd.bf/ntotal) * z.bft
        upper.bft <- pd + sqrt(vd.bf/ntotal) * z.bft
        p.bfti <- pt(t.bf, df.sw)
        p.bft <- min(2 - 2 * p.bfti, 2 * p.bfti)},
    "greater"={z.bft <- qt(1 - conflevel, df.sw)
        lower.bft <- pd - sqrt(vd.bf/ntotal) * z.bft
        upper.bft <- 1
        p.bft <- 1 - pt(t.bf, df.sw)},
    "less"={z.bft <- qt(1 - conflevel, df.sw)
        upper.bft <- pd + sqrt(vd.bf/ntotal) * z.bft
        lower.bft <- 0
        p.bft <- pt(t.bf, df.sw)})


parameter <- df.sw
names(parameter)<-"df"
conf.int <- c(lower.bft, upper.bft)
t.value <- t.bf
names(t.value)<-"t.bf"
p.value <- p.bft
method <- "Relative Effects, t-Approximation"
    })

names(estimate) <- cmpid
attr(conf.int, which="conf.level") <- conf.level
null.value <- 0.5; names(null.value) <- "relative effect"

out<-list(
statistic=t.value,
parameter=parameter,
p.value=p.value,
conf.int=conf.int,
estimate=estimate,
null.value=null.value,
alternative=alternative,
method=method,
data.name=dname)

class(out)<-"htest"
return(out)

}