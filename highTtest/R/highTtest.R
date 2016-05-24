highTtest <- function(dataSet1, 
                      dataSet2,  
                      gammas,  
                      compare="BOTH",
                      cSequence = NULL,
                      tSequence = NULL){

  compare <- toupper(compare)

  if("BOTH" %in% compare) compare <- c("BH","ST")

  if(!all(compare %in% c("BH","ST","NONE"))) {
    stop("method specified is not available")
  }

  if(is.data.frame(dataSet1)){
    dataSet1 <- data.matrix(dataSet1)
  } else if(!is.matrix(dataSet1)){
    dataSet1 <- matrix(dataSet1,nrow=1)
  }

  if(is.data.frame(dataSet2)){
    dataSet2 <- data.matrix(dataSet2)
  } else if(!is.matrix(dataSet2)){
    dataSet2 <- matrix(dataSet2,nrow=1)
  }

  if(nrow(dataSet2) != nrow(dataSet1)){ 
    stop("data sets must have same # of observations")
  }

  n1 <- ncol(dataSet1)
  n2 <- ncol(dataSet2)

  Tstar <- tStatistic(data1 = dataSet1,
                      data2 = dataSet2)

  p_value <- Tstar$pv
  po <- order(p_value)
  p_value <- p_value[po]
  rm(dataSet1, dataSet2)

  tmp <- CaoKosorok(Tstar = Tstar$statistic, 
                    gammas = gammas,
                    cSequence = cSequence,
                    tSequence = tSequence)
  indicator.ck <- tmp$ind
  pi1 <- tmp$pi1

  if("BH" %in% compare){
    indicator.bh <- BenjaminiHochberg(alphas = gammas, p_value = p_value)
    indicator.bh[po,] <- indicator.bh
  } else {
    indicator.bh <- NULL
  }

  if("ST" %in% compare){
    indicator.st <- StoreyTibshirani(alphas = gammas, p_value = p_value)
    indicator.st[po,] <- indicator.st
  } else {
    indicator.st <- NULL
  }

  p_value[po] <- p_value

  HCO <- new("highTtest",
             CK = indicator.ck,
             pi1 = pi1,
             pvalue = p_value,
             BH = indicator.bh,
             ST = indicator.st,
             gammas = gammas)


  return(HCO)

}

