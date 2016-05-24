#función matriz.mpa, a partir del archivo leído mediante leer.mpa, crea la matriz
# de coeficientes de asociación E, y una lista de las palabras
#parámetros: leer.mpa: el objeto que resulta de la funcion leer.mpa
#            sep.ind: el separador de individuos
#            sep.pal: el separador de palabras
#            fmin: frecuencia mínima de las palabras
#            cmin: co-ocurrecia mínima entre dos palabras

matriz.mpa<-function(leer.mpa, sep.ind="ind0", sep.pal="/", fmin=3, cmin=3)
{
  ocurrencias <- unlist(strsplit(leer.mpa[1:length(leer.mpa)],sep.pal))
  ocurrencias <- ocurrencias[nchar(ocurrencias) > 0]
  ocurrencias <- tolower(ocurrencias)

  linea <- rep(0, length.out = length(ocurrencias))
  which.index <- grep(sep.ind, ocurrencias)
  full.index <- which.index
     if (length(which.index) > 0) {
        which.index.2 <- c((which.index[2:length(which.index)] - 1), (length(ocurrencias)))
        for (i in 1:length(which.index)) {
            linea[which.index[i]:which.index.2[i]] <- i
        }
    }
  ocurrencias <- ocurrencias[-full.index]
  linea <- linea[-full.index]
  tl <- table(linea,ocurrencias)
  E1 <- t(tl)%*%tl
  E1 <- E1[colSums(tl)>=fmin,colSums(tl)>=fmin]
  E2 <- ifelse(E1<cmin,0,E1)
  di <- diag(E2)
  dm <- solve(diag(di))
  E <- dm%*%(E2*E2)%*%dm
  rownames(E) <- colnames(tl)[colSums(tl)>=fmin]
  colnames(E) <- colnames(tl)[colSums(tl)>=fmin]
  rownames(E1) <- colnames(tl)[colSums(tl)>=fmin]
  colnames(E1) <- colnames(tl)[colSums(tl)>=fmin]
  nombres <- colnames(tl)[colSums(tl)>=fmin]
  res <- list(Matriza = E, Matrizc = E1, Palabras = nombres, tl=tl)
  return(res)
}