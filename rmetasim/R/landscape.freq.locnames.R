landscape.freq.locnames <- function(Rland)
  {
      l <- Rland
      num.loc <- length(landscape.ploidy(l))
      namevec <- NULL
      for (loc in 1:num.loc)
      {
          genos <- landscape.locus(l,loc)[,-1:-landscape.democol()]
          loc.names <- paste(loc, names(table(unlist(genos))), sep = ".")
          namevec <- c(namevec,paste("L", loc.names, sep = '')) #not a fast construct, I know.  But remember Knuth "early optimization is the root of all
      }
      namevec
  }
