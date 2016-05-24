move <- function(from, to){
  noms <- as.list(match.call())
  if(noms$from != noms$to){
    eval.parent(parse(text = paste(noms$from, "->", noms$to)))
    eval.parent(parse(text = paste("rm(", noms$from, ")")))
  }
}

mv <- move
