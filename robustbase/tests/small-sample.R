library(robustbase)
## testing functions:
source(system.file("xtraR/ex-funs.R", package = "robustbase"))

set.seed(152)

Nmax <- 12
nn <- length(nset <- c(2:Nmax, 20, 50))## NOTA BENE: n == 1 etc are NOT YET TREATED!
Sim <- 2^9 # = 512

sn <- qn <- numeric(Sim)
cpu <- numeric(nn)
names(cpu) <- as.character(nset)

for(n in nset) {
    nS <- Sim ## if(n < 20) Sim else round(10*Sim/n)
    cat("\nn = ",n,"\n------\nno.Sim. = ",nS,"\n")
    cpu[as.character(n)] <- system.time(for(i in 1:nS) {
        x <- rnorm(n)
        sn[i] <- Sn0R(x)
        qn[i] <- Qn0R(x)
        Sn.x <- Sn(x, const = 1)
        Qn.x <- Qn(x, const = 1)
        if(!is.all.equal(Sn.x, sn[i], tol = 1e-5))
            cat("i=",i," Sn() != Sn0R():  ", Sn.x, "!=", sn[i],"\n")
        if(!is.all.equal(Qn.x, qn[i], tol = 1e-5))
            cat("i=",i," Qn() != Qn0R():  ", Qn.x, "!=", qn[i],"\n")
    })[1]
    cat("Mean and its. std.err;  Quartiles of Sn(x_1 .. x_n) and Qn(...):\n")
    print(c(mean(sn), sd(sn)/sqrt(nS), quantile(sn, p = (1:3)/4)))
    print(c(mean(qn), sd(qn)/sqrt(nS), quantile(qn, p = (1:3)/4)))
}

rbind("Time (CPU) used:" = summary(cpu))
