#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# plugInValue : Estimate plug-in value                                         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# optTx1   : vector of first-stage treatments corresponding to the first-stage #
#            decision rule of the proposed regime.                             #
# optTx2   : vector of second-stage treatments corresponding to the            #
#            second-stage decision rule of the proposed regime.                #
# response : vector of response                                                #
# tx1      : vector of first-stage randomized treatments                       #
# tx2      : vector of second-stage randomized treatments                      #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
plugInValue <- function(optTx1, 
                        optTx2, 
                        response, 
                        tx1, 
                        tx2){

  if( is(optTx1,"factor") ) optTx1 <- as.numeric(levels(optTx1)[optTx1])
  if( is(optTx2,"factor") ) optTx2 <- as.numeric(levels(optTx2)[optTx2])
  #--------------------------------------------------------------------------#
  # determine the average value for patients that followed the optimal regime#
  #--------------------------------------------------------------------------#
  tx1 <- as.integer(round(tx1,0L))
  tx2 <- as.integer(round(tx2,0L))
  optTx1 <- as.integer(round(optTx1,0L))
  optTx2 <- as.integer(round(optTx2,0L))

  ind <- as.numeric( tx1 == optTx1 ) * as.numeric( tx2 == optTx2 )
  value <- sum(response*ind)/sum(ind)

  txOptions1 <- sort(unique(c(tx1,optTx1)))
  txOptions2 <- sort(unique(c(tx2,optTx2)))

  #--------------------------------------------------------------------------#
  # determine average value for patients that followed all get 1             #
  #--------------------------------------------------------------------------#
  fixedReg <- matrix(data = 0,
                     nrow = length(txOptions1),
                     ncol = length(txOptions2),
                     dimnames = list(paste("tx1=", txOptions1, sep=""),
                                     paste("tx2=", txOptions2, sep="")))

  for( i in txOptions1 ) {
    for( j in txOptions2 ) {
      ind <- as.numeric(tx1 == i) * as.numeric(tx2 == j)
      fixedReg[i,j] <- sum(response * ind)/sum(ind)
    }
  }

  return(list(   "value" = value, 
              "fixedReg" = fixedReg))

}


