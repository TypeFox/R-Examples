makeDomEpi <- function(pedigree, output = c("AD", "DD", "both"), parallel = FALSE, invertD = FALSE, det = FALSE)
{
  type <- match.arg(output)
  Dout <- makeD(pedigree, parallel = parallel, invertD = invertD, returnA = TRUE)

  if(type == "AD"){
    AD <- Dout$A * Dout$D
    if(det) logDetAD <- determinant(AD, logarithm = TRUE)$modulus[1] else logDetAD <- NULL
    ADinv <- as(solve(AD), "dgCMatrix")
    listADinv <-sm2list(ADinv, rownames=pedigree[,1], colnames=c("row", "column", "ADinverse"))
    DD <- NULL
    logDetDD <- NULL
    DDinv <- NULL
    listDDinv <- NULL
    }      

  if(type == "DD"){
    DD <- Dout$D * Dout$D
    if(det) logDetDD <- determinant(DD, logarithm = TRUE)$modulus[1] else logDetDD <- NULL
    DDinv <- as(solve(DD), "dgCMatrix")
    listDDinv<-sm2list(DDinv, rownames=pedigree[,1], colnames=c("row", "column", "DDinverse"))
    AD <- NULL
    logDetAD <- NULL
    ADinv <- NULL
    listADinv <- NULL
   }

  if(type == "both"){
    AD <- Dout$A * Dout$D
    ADinv <- Matrix(solve(AD), sparse=TRUE)
    listADinv <-sm2list(ADinv, rownames=pedigree[,1], colnames=c("row", "column", "ADinverse"))
    ADinv <- as(ADinv, "dgCMatrix")
    DD <- Dout$D * Dout$D
    if(det){
      logDetAD <- determinant(AD, logarithm = TRUE)$modulus[1]
      logDetDD <- determinant(DD, logarithm = TRUE)$modulus[1]
    } else{ logDetAD <- logDetDD <- NULL}
    DDinv <- as(solve(DD), "dgCMatrix")
    listDDinv<-sm2list(DDinv, rownames=pedigree[,1], colnames=c("row", "column", "DDinverse"))
     }
return(list(D=Dout$D, logDetD = Dout$logDet, AD=AD, logDetAD = logDetAD, DD=DD, logDetDD = logDetDD, Dinv=Dout$Dinv, ADinv=ADinv, DDinv=DDinv, listDinv=Dout$listDinv, listADinv=listADinv, listDDinv=listDDinv))

}

