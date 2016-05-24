Addlist <- function(tree, string, word, node = 1) {
  b <- tree
  key <- str_sub(string,1,1)
  if (is.na(b$ch[node])) {
    b$ch[node] <- key
    b$flag[node] <- 0
  }
  if (key < b$ch[node]) {
    if (is.na(b$L[node])) {
      b$L[node] <- length(b$ch) + 1
      b <- Addlist(b, string, word, b$L[node])
    }
    else{
      b <- Addlist(b, string, word, b$L[node])
    }
  }
  else if (key > b$ch[node]) {
    if (is.na(b$R[node])) {
      b$R[node] <- length(b$ch) + 1
      b <- Addlist(b, string, word, b$R[node])
    }
    else{
      b <- Addlist(b, string, word, b$R[node])
    }
  }
  else {
    if (nchar(string) == 1) {
      b$flag[node] <- 1
      if(length(b$word) < node){b$word[[node]] <- word}
      else{b$word[[node]] <- append(b$word[[node]], word)}
    }
    else {
      if (!is.na(b$C[node]) || !is.na(b$L[node]) || !is.na(b$R[node])) {
        if (is.na(b$C[node]) && (b$flag == 1)) {
          b$C[node] <- length(b$ch) + 1
          b <- Addlist(b, str_sub(string,2,-1), word, b$C[node])
        }
        else if (is.na(b$C[node]) && (b$flag == 0)) {
          b$C[node] <- node + 1
          b <- Addlist(b, str_sub(string,2,-1), word, b$C[node])
        }
        else{
          b <- Addlist(b, str_sub(string,2,-1), word, b$C[node])
        }
      }
      else{
        b$C[node] <- length(b$ch) + 1
        b <- Addlist(b, str_sub(string,2,-1), word, b$C[node])
      }
    }
  }
  return(b)
}
