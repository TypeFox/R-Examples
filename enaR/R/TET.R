#' TET.R  --- TOTAL ENVIRON THROUGHFLOW
#' INPUT = network model
#' OUTPUT = total environ throughput - unit and scaled
#'
#' Borrett | July 7, 2012
#' ---------------------------------------------------

TET <- function(x,balance.override=FALSE){

                                        #Check for network class
  if (class(x) != 'network'){warning('x is not a network class object')}

                                        #Check for balancing
  if (balance.override){}else{
    if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' = ssCheck(x)}
    if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
  }

  oo <- get.orient() #original orientation
  if (oo == 'school'){oo <- 'internal'}
  set.orient('internal')
  E <- enaEnviron(x)
  set.orient(oo)
  input <- unpack(x)$z   # get data input elements
  output <- unpack(x)$y  # get data output elements

  # UNIT
  unit.input <- 0
  unit.output <- 0

  # REALIZED
  realized.input <- 0 # initialize
  realized.output <- 0 # initialize

  # UNIT & SCALED
  for(i in 1:length(input)){
    realized.input[i] = -sum(diag(E$input[[i]] * output[i]))
    realized.output[i] = -sum(diag(E$output[[i]] * input[i]))
    unit.input[i] = -sum(diag(E$input[[i]]))
    unit.output[i] = -sum(diag(E$output[[i]]))
  }

  return(
         list("realized.input"=realized.input,
              "realized.output"=realized.output,
              "unit.input"=unit.input,
              "unit.output"=unit.output)
         )

}
