is.ambig <-
function(DNAbin){
   x <- as.matrix(DNAbin)
   bases <- c(136, 72, 40, 24)
   ambig <- apply(x, 2, FUN=function(x) sum(as.numeric(!as.numeric(x) %in% bases)))
   ambig > 0
}

