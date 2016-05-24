spdfs <- function(tree,nodo,match) {
  b<-tree
  guardo2<-c()
  if (!is.na(nodo) && b$flag[nodo]==1) {
    guardo2 <- append(guardo2, match)
  }
  if (is.na(b$C[nodo]) && is.na(b$L[nodo]) && is.na(b$R[nodo])) {
  }
  if (!is.na(b$C[nodo])) {
    guardo2<- append(guardo2,spdfs(b,b$C[nodo],paste0(match,b$ch[b$C[nodo]])))
  }
  if (!is.na(b$R[nodo])) {
    guardo2<-append(guardo2,spdfs(b,b$R[nodo],paste0(str_sub(match,1,-2),b$ch[b$R[nodo]])))
  }
  if (!is.na(b$L[nodo])) {
    guardo2<-append(guardo2,spdfs(b,b$L[nodo],paste0(str_sub(match,1,-2),b$ch[b$L[nodo]])))
  }
  return(guardo2)
}