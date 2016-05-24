library(NPMVCP)
data("Alsmelterdata")
output <- NPMVCP(Alsmelterdata)

# p = dimension of each data vector #
p <- ncol(Alsmelterdata)

# c = degree of quarantine #
c <- 15

# N = total number of observation vectors #
N <- nrow(Alsmelterdata)

# set monitoring start value #
monitoring.start <- max(p + 10, 2*c + 3)

# load control limits #
CLdatastring <- paste("CLp", p, "c", c, sep="")
data(list=CLdatastring)
CL <- get(CLdatastring)[,"0.002"]

# extrapolate control limits beyond n = 500, if necessary #
nmax <- N-(monitoring.start-1)
if (nmax > 500) {
    ninv <- 1/(100:500)
    CLexmodel <- lm(CL[100:500] ~ ninv)
    CL <- c(CL, pmax(CL[500], CLexmodel$coef[1] + CLexmodel$coef[2]*(1/(501:nmax))))
}

# adjust index of control limit vector for plotting #
CLtoplot <- c(rep(NA, monitoring.start-1), CL[1:(N-(monitoring.start-1))])

\dontrun{plot(output$Rmax, ylim=c(0, 24), type="b")}
\dontrun{lines(CLtoplot)}
