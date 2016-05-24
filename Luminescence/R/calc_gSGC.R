#' Calculate De value based on the gSGC by Li et al., 2015
#'
#' Function returns De value and De value error using the global standardised growth
#' curve (gSGC) assumption proposed by Li et al., 2015 for OSL dating of sedimentary quartz
#'
#' The error of the De value is determined using a Monte Carlo simulation approach.
#' Solving of the equation is realised using \code{\link{uniroot}}.
#' Large values for \code{n.MC} will significantly increase the computation time.
#'
#'
#' @param data \code{\link{data.frame}} (\bold{required}): input data of providing the following
#' columns: 'LnTn', 'LnTn.error', Lr1Tr1', 'Lr1Tr1.error', 'Dr1'
#' Note: column names are not required. The function expect the input data in the given order
#'
#' @param gSGC.type \code{\link{character}} (with default): define the function parameters that
#' should be used for the iteration procedure: Li et al., 2015 (Table 2)
#' presented function parameters for two dose ranges: \code{"0-450"} and \code{"0-250"}
#'
#' @param gSGC.parameters \code{\link{list}} (optional): option to provide own function
#' parameters used for #' fitting as named list.
#' Nomenclature follows Li et al., 2015, i.e.
#' \code{list(A,A.error,D0,D0.error,c,c.error,Y0,Y0.error,range)}, range requires a vector for
#' the range the function is considered as valid, e.g. \code{range = c(0,250)}\cr
#' Using this option overwrites the default parameter list of the gSGC, meaning the argument
#' \code{gSGC.type} will be without effect
#'
#' @param n.MC \code{\link{integer}} (with default): number of Monte Carlo simulation runs for
#' error estimation, s. details.
#'
#' @param verbose \code{\link{logical}}: enable or disable terminal output
#'
#' @param plot \code{\link{logical}}: enable or disable graphical feedback as plot
#'
#' @param ... parameters will be passed to the plot output
#'
#' @return Returns an S4 object of type \code{\linkS4class{RLum.Results}}.
#' Slot \code{data} contains a \code{\link{list}} with the following structure:\cr
#' $ De.value (data.frame) \cr
#'  .. $ De  \cr
#'  .. $ De.error \cr
#'  .. $ Eta \cr
#' $ De.MC (list) contains the matricies from the error estimation.\cr
#' $ uniroot (list) contains the uniroot outputs of the De estimations
#' $ call (call) the original function call
#'
#'
#' @section Function version: 0.1.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montagine (France)\cr
#'
#' @seealso \code{\linkS4class{RLum.Results}}, \code{\link{get_RLum}}, \code{\link{uniroot}}
#'
#' @references  Li, B., Roberts, R.G., Jacobs, Z., Li, S.-H., 2015. Potential of establishing
#' a 'global standardised growth curve' (gSGC) for optical dating of quartz from sediments.
#' Quaternary Geochronology 27, 94-104. doi:10.1016/j.quageo.2015.02.011
#'
#' @keywords datagen
#'
#' @examples
#' results <- calc_gSGC(data = data.frame(
#' LnTn =  2.361, LnTn.error = 0.087,
#' Lr1Tr1 = 2.744, Lr1Tr1.error = 0.091,
#' Dr1 = 34.4))
#'
#' get_RLum(results, data.object = "De")
#'
#' @export
calc_gSGC<- function(
  data,
  gSGC.type = "0-250",
  gSGC.parameters,
  n.MC = 100,
  verbose = TRUE,
  plot = TRUE,
  ...
){

##============================================================================##
##CHECK INPUT DATA
##============================================================================##

  if(!is(data, "data.frame")){stop("'data' needs to be of type data.frame.")}
  if(!is(gSGC.type, "character")){stop("'gSGC.type' needs to be of type character.")}

  ##check length of input data
  if(ncol(data) != 5){stop("Structure of 'data' does not fit the expectations.")}

  ##rename columns for consistency reasons
  colnames(data) <- c('LnTn', 'LnTn.error', 'Lr1Tr1', 'Lr1Tr1.error', 'Dr1')


##============================================================================##
##DEFINE FUNCTION
##============================================================================##

    ##define function, nomenclature according to publication that should be solved
    f <- function(x,A,D0,c,Y0,Dr1,Lr1Tr1,LnTn) {
      (((A * (1 - exp( - Dr1 / D0))) + c * Dr1 + Y0)/Lr1Tr1) -
      (((A * (1 - exp( - x/D0))) + c * x + Y0)/LnTn)
  }

    ##set general parameters
    if (!missing(gSGC.parameters)) {
      A <- gSGC.parameters$A
      A.error <- gSGC.parameters$A.error
      D0 <- gSGC.parameters$D0
      D0.error <- gSGC.parameters$D0.error
      c <- gSGC.parameters$c
      c.error <- gSGC.parameters$c.error
      Y0 <- gSGC.parameters$Y0
      Y0.error <- gSGC.parameters$Y0.error
      range <- gSGC.parameters$range

    }else{
      if (gSGC.type == "0-450") {
        A <- 0.723
        A.error <- 0.014
        D0 <- 65.1
        D0.error <- 0.9
        c <- 0.001784
        c.error <- 0.000016
        Y0 <- 0.009159
        Y0.error <- 0.004795

        range <- c(0.1,250)

      }else if (gSGC.type == "0-250") {
        A <- 0.787
        A.error <- 0.051
        D0 <- 73.9
        D0.error <- 2.2
        c <- 0.001539
        c.error <- 0.000068
        Y0 <- 0.01791
        Y0.error <- 0.00490

        range <- c(0.1,250)

      }else{
        stop("Unknown input for 'gSGC.type'")

      }

    }

    ##Define size of output objects
    output.data <- data.table(
      De = numeric(length = nrow(data)),
      De.error =  numeric(length = nrow(data)),
      Eta =  numeric(length = nrow(data))
    )

    ##set list for De.MC
    output.De.MC <- vector("list",  nrow(data))

    ##set list for uniroot
    output.uniroot <-  vector("list",  nrow(data))


##============================================================================##
##CALCULATION
##============================================================================##


 for(i in 1:nrow(data)){

    Lr1Tr1 <-data[i,"Lr1Tr1"]
    Lr1Tr1.error <- data[i,"Lr1Tr1.error"]
    Dr1 <- data[i,"Dr1"]
    Dr1.error <- data[i,"Dr1.error"]

    LnTn <- data[i,"LnTn"]
    LnTn.error <- data[i,"LnTn.error"]

  ##calculate mean value
    temp <- try(uniroot(
      f,
      interval = c(0.1,450),
      tol = 0.001,
      A = A,
      D0 = D0,
      c = c,
      Y0 = Y0,
      Dr1 = Dr1,
      Lr1Tr1 = Lr1Tr1,
      LnTn = LnTn,
      extendInt = 'yes',
      check.conv = TRUE,
      maxiter = 1000
    ), silent = TRUE)

  if(!inherits(temp, "try-error")){

    ##get De
    De <- temp$root

    ##calculate Eta, which is the normalisation factor
    Eta <- ((A * (1 - exp( - Dr1 / D0))) + c * Dr1 + Y0)/Lr1Tr1

    ##--------------------------------------------------------------------------##
    ##Monte Carlo simulation for error estimation

    ##set matrix
    temp.MC.matrix <- matrix(nrow = n.MC, ncol = 8)

    ##fill matrix
    temp.MC.matrix[,1:6] <- matrix(rnorm(
      n.MC * 6,
      mean = c(LnTn, Lr1Tr1, A, D0, c, Y0),
      sd = c(LnTn.error, Lr1Tr1.error, A.error, D0.error, c.error, Y0.error)
    ), ncol = 6, byrow = TRUE)


      ##run uniroot to get the De
      temp.MC.matrix[,7] <- sapply(1:n.MC, function(x){

        uniroot(f,
                interval = c(0.1,450),
                tol = 0.001,
                A = temp.MC.matrix[x,3],
                D0 = temp.MC.matrix[x,4],
                c = temp.MC.matrix[x,5],
                Y0 = temp.MC.matrix[x,6],
                Dr1 = Dr1,
                Lr1Tr1 =temp.MC.matrix[x,2],
                LnTn = temp.MC.matrix[x,1],
                check.conv = TRUE,
                extendInt = 'yes',
                maxiter = 1000
                )$root

      })

      ##calculate also the normalisation factor
      temp.MC.matrix[,8] <- (temp.MC.matrix[,3] * (1 - exp( - Dr1 / temp.MC.matrix[,4])) +
        temp.MC.matrix[,5] * Dr1 + temp.MC.matrix[,6])/temp.MC.matrix[,2]


      ##re-name matrix
      colnames(temp.MC.matrix) <- c("LnTn","Lr1Tr1","A","D0","c","Y0","De","Eta")

      ##get De error as SD
      De.error <- sd(temp.MC.matrix[,7])

  }else{
    warning("No solution was found!")
    De <- NA
    Eta <- NA
    De.error <- NA

    ##set matrix
    temp.MC.matrix <- matrix(nrow = n.MC, ncol = 8)

    ##fill matrix
    temp.MC.matrix[,1:6] <- matrix(rnorm(
      n.MC * 6,
      mean = c(LnTn, Lr1Tr1, A, D0, c, Y0),
      sd = c(LnTn.error, Lr1Tr1.error, A.error, D0.error, c.error, Y0.error)
    ), ncol = 6, byrow = TRUE)


  }

##============================================================================##
##PLOT OUTPUT
##============================================================================##

  if (plot) {

    ##set plot settings
    plot.settings <- list(
      main = "gSGC and resulting De",
      xlab = "Dose/(a.u.)",
      ylab = expression(paste("Re-norm. ", L[x]/T[x])),
      xlim = NULL,
      ylim = NULL,
      lwd = 1,
      lty = 1,
      pch = 1,
      col = "red",
      grid = expression(nx = 10, ny = 10),
      mtext = ""
    )

    plot.settings <-  modifyList(plot.settings, list(...))



    ##graphical feedback
    x <- NA
    curve(
      A * (1 - exp(-x / D0)) + c * x + Y0, from = 0, to = 500,
      xlab = plot.settings$xlab,
      ylab = plot.settings$ylab,
      main = plot.settings$main,
      xlim = plot.settings$xlim,
      ylim = plot.settings$ylim,
      lwd = plot.settings$lwd,
      lty = plot.settings$lty
    )

    mtext(side = 3, plot.settings$mtext)

    if(!is.null(plot.settings$grid)){
      graphics::grid(eval(plot.settings$grid))

    }

    points(temp$root,Eta*LnTn, col = plot.settings$col, pch = plot.settings$pch)
    segments(De - De.error,Eta * LnTn,
             De + De.error,Eta * LnTn)

    hist(temp.MC.matrix[,7], freq = FALSE,add = TRUE)


  }

##============================================================================##
##OUTPUT VISUALISATION
##============================================================================##

    if (verbose) {
      cat("\n[calc_gSGC()]")
      cat("\n\t Corresponding De based on the gSGC\n")

      cat(paste0("\n\t"," Ln/Tn:\t\t ",LnTn," \u00B1 ", LnTn.error,"\n"))
      cat(paste0("\t"," Lr1/Tr1:\t ",Lr1Tr1," \u00B1 ", Lr1Tr1.error,"\n"))
      cat(paste0("\t"," Dr1:\t\t ",Dr1,"\n"))
      cat(paste0("\t"," f(D):\t\t ",A," * (1 - exp(-D /",D0,")) + c * D + ",Y0,"\n"))
      cat(paste0("\t"," n.MC:\t\t ",n.MC,"\n"))
      cat(paste0("\t ------------------------------ \n"))
      cat(paste0("\t De:\t\t",round(De,digits = 2)," \u00B1 ",round(De.error,digits = 2),"\n"))
      cat(paste0("\t ------------------------------ \n"))

    }


##============================================================================##
##CREATE OUTPUT OBJECTS
##============================================================================##

    ##needed for data.table
    temp.De <- De
    temp.De.error <- De.error
    temp.Eta <- Eta

    ##replace values in the data.table with values
    output.data[i, `:=` (De = temp.De,
                         De.error = temp.De.error,
                         Eta = temp.Eta)]

    rm(list = c('temp.De', 'temp.De.error', 'temp.Eta'))

    ##matrix - to prevent memory overload limit output
    if(n.MC * nrow(data) > 1e6){

      if(i == 1){

        output.De.MC[[i]] <- temp.MC.matrix

      }else{

        output.De.MC[[i]] <- NA

      }

      warning("Only the first MC matrix is returned to prevent memory overload!")

    }else{

      output.De.MC[[i]] <- temp.MC.matrix

    }


    output.uniroot[[i]] <- temp


}##end for loop

##============================================================================##
##OUTPUT RLUM
##============================================================================##

    temp.RLum.Results <- set_RLum(
      class = "RLum.Results",
      data = list(
        De = as.data.frame(output.data),
        De.MC =  output.De.MC,
        uniroot = output.uniroot,
        call = sys.call()
      )
    )

  return(temp.RLum.Results)
}
