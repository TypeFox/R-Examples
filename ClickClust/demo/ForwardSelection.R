
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

set.seed(1234)
      
# FORWARD STATE SELECTION
     
F2 <- click.forward(X = C$X, K = 2)
print.search(F2)
