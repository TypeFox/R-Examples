# -------------------------------------------
# Funcion para calcular centroides de una particion
# codigo tomado de funcion FactoClasss.class de ade4
# julio 7 de 2009
# Campo E. pardo
# Entra particion, coordenadas y pesos
# salen coordenadas de los centroides
# ---------------------------------------------
centroids <- function(df,cl,rw=rep(1/nrow(df),nrow(df)))
{
  f1 <- function(cl)
  {
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    data.frame(x)
  }
  f1(cl)  # TDC de la particion
  dfdistri <- f1(cl) * rw  
  w1 <- unlist(lapply(dfdistri, sum)) # pesos de las clases
  coo <- as.matrix(t(dfdistri)) %*%as.matrix(df)/w1
    if (class(cl)=="factor") rownames(coo) <-levels(cl)
    return(list(centroids=coo,weights=w1))
}                                         


