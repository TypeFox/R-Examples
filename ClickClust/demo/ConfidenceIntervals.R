
data(synth)
     
# FUNCTION THAT REPLACES CHARACTER STATES WITH NUMERIC VALUES
repl.levs <- function(x, ch.lev){
        for (j in 1:length(ch.lev)) x <- gsub(ch.levs[j], j, x)
        return(x)
}
     
# DETECT ALL STATES IN THE DATASET
d <- paste(synth$data, collapse = " ")
d <- strsplit(d, " ")[[1]]
ch.levs <- levels(as.factor(d))
     
# CONVERT DATA TO THE FORM USED BY click.read()
S <- strsplit(synth$data, " ")
S <- sapply(S, repl.levs, ch.levs)
S <- sapply(S, as.numeric)

# READ DATA
C <- click.read(S)

set.seed(123)

# EM ALGORITHM (without initial state probabilities)
     
N2 <- click.EM(X = C$X, K = 2)

V <- click.var(X = C$X, alpha = N2$alpha, gamma = N2$gamma, z = N2$z)
st.err <- sqrt(diag(V))

Estimates <- c(N2$alpha[-2], as.vector(apply(N2$gamma[,-5,], 3, t)))

# 95% confidence intervals for parameter estimates

Lower <- Estimates - qnorm(0.975) * st.err
Upper <- Estimates + qnorm(0.975) * st.err

cbind(Estimates, Lower, Upper)
