
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

# EM ALGORITHM (with initial state probabilities)
     
M2 <- click.EM(X = C$X, y = C$y, K = 2)

click.plot(X = C$X, y = C$y, id = M2$id, states = ch.levs, obs.lwd = 0.3)
