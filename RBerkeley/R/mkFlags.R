mkFlags <- function(...) {
  if(nargs()==1) {
    flags <- gsub(" ","",
                  unlist(
                    strsplit(
                      as.character(
                        match.call(call=sys.call())[-1]),"\\|")))
    flags <- gsub('\"',"",flags)
  } else {
    flags <- as.character(match.call(call=sys.call())[-1])
  }
  .Call("mkFlags", flags)
}
