## Loop through all attached directories looking for
## regular expression pattern.

objip <-
  function(pattern, where = search(), all.names=FALSE, mode="any", class,  ## , sorted=TRUE (requires R-3.2.0)a
           ls.function=if (mode != "any" || !missing(class)) "ls.str" else "ls")
    {
      ls.function <- match.arg(ls.function, c("ls", "ls.str"))
      result <- list()
      for(i in match(where, search())) {
        obj <-
          if (ls.function=="ls")
            ls(pos=i, pattern = pattern, all.names=all.names) ## , sorted=sorted
          else
            ls.str(pos=i, pattern = pattern, all.names=all.names, mode=mode)
        if(length(obj) > 0)
          result[[where[i]]] <- obj
      }
      if (ls.function=="ls.str" && !missing(class))
        for (i in names(result)) {
          keep <- sapply(result[[i]], function(x, class) is(get(x), class), class)
          result[[i]] <- result[[i]][keep]
          if (length(result[[i]]) == 0) result[[i]] <- NULL
        }
      result
    }


if (FALSE) {
  objip(pat="AE")
  objip(pat="AE",  class="data.frame")
  objip(pat="AE",  mode="function")
  objip(pat="AE",  class="function")
}
