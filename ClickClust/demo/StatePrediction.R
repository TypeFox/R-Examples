
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

T <- click.predict(M = 3, gamma = N2$gamma, pr = N2$z[1,])
colnames(T) <- ch.levs
T[S[[1]][length(S[[1]])],]

