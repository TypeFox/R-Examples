crearFullTree <- function(input,maxdist){
  if(maxdist==1){
    arbol<-spellTree_1()
  }
  else if(maxdist==2){
    arbol<-spellTree_2()
  }
  else if(maxdist==3){
    arbol<-spellTree_3()
  }
  else{stop("argument \"maxdist\" out of range (1:3)")}
  contador <- 0
  if(!file.exists(input)[1]){
    lista<-delList(input,maxdist)
    for(i in 1:length(lista)){
      word<-names(lista[i])
      for (item in lista[[i]]) {
        arbol <- Addlist(arbol, as.character(item),word)
        contador <- contador + 1
        if (contador %% 10000 == 0) {
          print(contador)
        }
      }
    }
  }
  else if (is.null(input)) {return()}
  else{
    content <- readLines(input, warn = FALSE)
    lista<-delList(content,maxdist)
    for(i in 1:length(lista)){
      word<-names(lista[i])
      for (item in lista[[i]]) {
        arbol <- Addlist(arbol, as.character(item),word)
        contador <- contador + 1
        if (contador %% 10000 == 0) {
          print(contador)
        }
      }
    }
  }
#   cat(paste("Tree created with",contador,"words and",length(arbol$ch),"nodes.\n"))
  return(arbol)
}
