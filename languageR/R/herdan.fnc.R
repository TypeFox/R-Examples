`herdan.fnc` <-
function(text, chunks) {
  nchunks = length(chunks)
  tokens = chunks
  types = rep(0, nchunks)
  if (length(chunks) < 2) stop("number of chunks should be at least 2")
  if (length(chunks) != length(unique(chunks))) 
    stop("duplicate values in chunks")
  for (i in 1:nchunks) {
    types[i] = length(unique(text[1:chunks[i]]))
  }
  return(list(growth = data.frame(tokens, types),
              C = coef(lm(log(types)~log(tokens)))[2]))
}

