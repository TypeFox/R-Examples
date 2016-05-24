tableToDescription2 <-
function(data) {
  data[,1:3] -> model
  as.data.frame(data[,4]) -> spp.data
  spaces <-function(x) return(gsub("(?<=[\\s])\\s*|^\\s+$", "", x, perl=TRUE))
    for (k in 1:nrow(data)) {
      model[k,] -> m0
      as.matrix(m0) -> m0
      m0[,3] -> sep0
      m0[,1] -> pre0
      m0[,2] -> suf0
      if (pre0 == "") {
        pre0[] <- "xxxwww"
      }
      if (suf0 == "") {
        suf0[] <- "xxxwww"
      }
      as.character(spp.data[k,1]) -> d0
      paste(pre0, " ", d0, " ", suf0, sep0, " ", sep="") -> s0
      sub("xxxwww ", "", s0) -> s0
      sub("xxxwww ", "", s0) -> s0
      sub(" xxxwww", "", s0) -> s0
      spaces(s0) -> s0
      cat(s0, file="", fill=F, append=T)
    }
  cat(" ", fill=T)
}
