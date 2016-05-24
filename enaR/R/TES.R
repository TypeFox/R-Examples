#' TES.R  --- TOTAL ENVIRON STORAGE
#' INPUT = network model
#' OUTPUT = total environ throughput - unit and scaled
#'
#' Borrett | July 7, 2012
#' ---------------------------------------------------

TES <- function(x,balance.override=FALSE){

                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' = ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }
                                        #
  oo <- get.orient() #original orientation
  if (oo == 'school'){oo <- 'internal'}
  set.orient('internal')
  S <- enaStorage(x)
  set.orient(oo)
  input <- unpack(x)$z   # get data input elements
  output <- unpack(x)$y  # get data output elements

  # UNIT
  X = S$S 
  unit.output <- apply(X,2,sum)
  X =  S$SP
  unit.input <- apply(X,2,sum)

  # REALIZED
  X = S$S %*% diag(input)
  realized.output <- apply(X,2,sum)
  X =  diag(output) %*% S$SP
  realized.input <- apply(X,2,sum)
  
  return(list("realized.input"=realized.input,"realized.output"=realized.output,"unit.input"=unit.input,"unit.output"=unit.output))
}
