"lcomoms2" <-
function(DATAFRAME, nmom=3, asdiag=FALSE, opdiag=FALSE, ...) {
  if(length(names(DATAFRAME)) != 2) {
    warning("a data.frame having only two columns is required")
    return(NULL)
  }
  if(nmom < 1) {
    warning("nmom is < 1 must be [1,5]")
    return(NULL)	
  }
  if(nmom > 5) {
    warning("nmom is > 5 must be [1,5]")
    return(NULL)	
  }
  
  L1 <- L2 <- L3 <- L4 <- L5 <- list()
  L1$matrix <- L2$matrix <- L3$matrix <-
               L4$matrix <- L5$matrix <- NULL
  T2 <- T3 <- T4 <- T5 <- list()
  T2$matrix <- T3$matrix <- T4$matrix <- T5$matrix <- NULL
  
  if(nmom >= 1) L1 <- Lcomoment.matrix(DATAFRAME,k=1)
  if(nmom >= 2) L2 <- Lcomoment.matrix(DATAFRAME,k=2)
  if(nmom >= 3) L3 <- Lcomoment.matrix(DATAFRAME,k=3)
  if(nmom >= 4) L4 <- Lcomoment.matrix(DATAFRAME,k=4)
  if(nmom >= 5) L5 <- Lcomoment.matrix(DATAFRAME,k=5)
  
  if(nmom >= 2) T2 <- Lcomoment.coefficients(L2,L2)
  if(nmom >= 3) T3 <- Lcomoment.coefficients(L3,L2)
  if(nmom >= 4) T4 <- Lcomoment.coefficients(L4,L2)
  if(nmom >= 5) T5 <- Lcomoment.coefficients(L5,L2)

  if(opdiag) {
    z <- list(L1=c(L1$matrix[1,2],L1$matrix[2,1]),
              L2=c(L2$matrix[1,2],L2$matrix[2,1]),
              T2=c(T2$matrix[1,2],T2$matrix[2,1]),
              T3=c(T3$matrix[1,2],T3$matrix[2,1]),
              T4=c(T4$matrix[1,2],T4$matrix[2,1]),
              T5=c(T5$matrix[1,2],T5$matrix[2,1]))
  } else if(asdiag) {
    z <- list(L1=diag(L1$matrix),
              L2=diag(L2$matrix),
              T2=diag(T2$matrix),
              T3=diag(T3$matrix),
              T4=diag(T4$matrix),
              T5=diag(T5$matrix))
  }
  else {
    z <- list(L1=L1$matrix,
              L2=L2$matrix,
              T2=T2$matrix,
              T3=T3$matrix,
              T4=T4$matrix,
              T5=T5$matrix)  	
  }
  z$nmom <- nmom
  z$source <- "lcomoms2"
  return(z)
}
