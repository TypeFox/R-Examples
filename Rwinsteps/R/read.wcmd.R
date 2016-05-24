read.wcmd <- function(filename) {

  cmdr <- readLines(filename)

  cmdr <- trim(gsub(";.+", "", cmdr))

  asti <- grep("\\*", cmdr)
  nami <- grep("=.+", cmdr, ignore.case = TRUE)
  nami <- nami[!nami %in% asti]
  endi <- grep("&end", cmdr, ignore.case = TRUE) + 1

  cmdnames <- trim(tolower(gsub("=.+", "", cmdr[nami])))
  out <- list()
  out[cmdnames] <- trim(tolower(gsub(".+=", "", cmdr[nami])))

  if(length(asti) > 0) {
    astnami <- asti[seq(1, length(asti), 2)]
    astendi <- asti[seq(2, length(asti), 2)]
    astnames <- trim(tolower(gsub("=.+", "",
      cmdr[astnami])))
    for(i in 1:length(astnames))
      out[[astnames[i]]] <- trim(tolower(cmdr[(astnami[i] + 1):
        (astendi[i] - 1)]))
  }

  if(any(tolower(cmdr) == "end labels")){
    out$labels <- cmdr[endi:(length(cmdr) - 1)]
    if(length(out$labels) != out$ni) {
      warning("Number of items, ni = ", out$ni,
        ", does not equal length of labels,",
        length(out$labels))
    }
  }

  out <- as.wcmd(out)

  return(out)
}
