####################################################################################################################
AQSys.mathDesc <- function(mathDesc) {
  # Each switch option provides an equation that will be available to be used to
  # make plots and calculate compositions for a system with known parameters
  switch(
    mathDesc,
    "merchuk" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        # merchuk's equation
        P1 * exp(P2 * (XC ^ (0.5)) - P3 * (XC ^ 3))
      }
    },
    "murugesan" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        # murugesan's equation
        P1 + P2 * (XC) ^ 0.5 + P3 * XC
      }
    },
    "tello" = {
      Fn <- function(CoefSET, XC) {
        # equation's parameters
        P1 <- CoefSET[1]
        P2 <- CoefSET[2]
        P3 <- CoefSET[3]
        # tello's equation
        P1 * log(P2 + XC) + P3
      }
    },
    # if user selects an option not available, it triggers an error
    # (check AQSys.err.R for details)
    AQSys.err("0")
  )
  # return chosen function
  return(Fn)
}
####################################################################################################################
#'@rdname AQSysList
#'@export AQSysList
#'@title Aqueous Systems Descriptors already implemented
#'@description The function returns a list of all mathematical descriptors available at the time.
AQSysList <- function() {
  # a new entry in updte must be added for each new equation implemmented in AQSys.mathDesc()
  # updte entries' name must match AQSys.mathDesc switch options
  updte <- c("merchuk", "murugesan", "tello")
  # return list
  return(updte)
}
####################################################################################################################
