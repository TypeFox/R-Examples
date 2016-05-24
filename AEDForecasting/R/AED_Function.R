#' CPI Function
#'
#' Incorporate change point analysis in ARIMA forecasting
#' @param  myts a time series object
#' @param  startChangePoint  a positive integer for minimum number of changepoints
#' @param  endChangePoint  a positive integer for maximum number of change points. If 0 then  only startChangePoint number of change points will be entered. Should be either 0 or greater than startChangePoint and if so the algorithm will loop through all values inbetween subject to step
#' @param  step  an integer to step through loop of change points
#' @param  num Bump model number (see below)
#' @param  cpmeth changepoint method. Default is BinSeg. See cpa package for details
#' @param  CPpenalty default is SIC. See cpa package for details
#' @param  showModel default is False, if True shows all models for all changepoints, if an integer all models for that changepoint, if a string all changepoints for that model
#' @return A data frame with all the results from analysis
#' @export
#' @importFrom "grDevices" "dev.new"
#' @importFrom "graphics" "grid" "plot"
#' @importFrom "stats" "is.ts"

cpi <- function(myts, startChangePoint = 1, endChangePoint = 0, step = 1, num=15, cpmeth='BinSeg', CPpenalty="SIC", showModel=FALSE) {

  if (is.ts(myts) == FALSE) {
    message("Error: First parameter must be a time series")
  } else if (!((startChangePoint > 0) && (startChangePoint %% 1 == 0))) {
    message("Error: Starting change point must be a positive integer")
  } else if (!((endChangePoint > 0) && (endChangePoint %% 1 == 0)) && (endChangePoint != 0)){

    message("Error: Ending change point must be a positive integer")
  } else if ((endChangePoint <= startChangePoint) && (endChangePoint != 0)) {
    message("Error: Ending change point must be greater than starting change point")
  } else {
       f <- function(t) { return (exp(-1/t) * (t > 0)) } # match the post
       StackCutoff <- function(t, start_level=1, start_transition=1, stop_transition=2) {
            tr <- 1 + (t - start_transition) / (stop_transition - start_transition)
            f2t <- f(2 - abs(tr))
            ft1 <- f(abs(tr) - 1)
            return (start_level * f2t / (ft1 + f2t))
       }

       ExpCos <- function(t, start_point, amplitude=1, decay.speed=0.1, period=20 ) {
            return ( amplitude * cos(2*pi*(t-start_point)/period) / exp((t-start_point) * decay.speed) )
       }

       ModulatedExp <- function(t, start_point, amplitude=0.1, decay.speed=0.1, period=20) {
            return ( (1 + amplitude * cos(2*pi*(t-start_point)/period) ) / exp((t-start_point) * decay.speed) )
       }
    DF2 <- 0
    DFtemp <- c()

    cpiP2 <- function(myts, n=startChangePoint, cpmeth.=cpmeth, CPpenalty.=CPpenalty, num.=num, showModel.=showModel) {

      DF <- data.frame('Model'=rep("", num), 'No.ChangePoints'=rep(n, num),
                       'AIC' = rep(0,num),stringsAsFactors = FALSE, 'cpiv' = I(vector(mode="list", length=num)))


      nSamples <- length(myts)
      model11_vector <- signal::bartlett(nSamples)
      model12_vector <- abs(signal::flattopwin(nSamples))
      showThisModel <- FALSE
      if (class(showModel) == "logical") { showThisModel <- ( showModel == TRUE) }

      for (k in 1:num) {
        #for (abc in 1: 2) {
           #k = 8

        m.bin <- suppressWarnings(changepoint::cpt.mean(myts,penalty=CPpenalty,method=cpmeth,Q=n))
        aveTS=mean(myts)
        temp1 <- c()
        temp2 <- c()
        cpiv <- matrix(0,ncol = n, nrow = length(myts))
        if (class(showModel) == "character") { showThisModel <- ( (sprintf("model %d", k) == showModel) ) }

        for (i in 1:n) {
          if (class(showModel) == "numeric") { showThisModel <- (showModel == i) }

          for (j in 1:m.bin@cpts[i]) {
            cpiv[j,i] <- 0 #for the ith interv var fills cpiv with 0s until the ith cpa row.
          }

          temp1 <- length(myts) - m.bin@cpts[i]

          for (j in 1:temp1) {
            temp2 <- j + m.bin@cpts[i]
            if(k == 1) {
              cpiv[temp2,i] <- 1 - (1-(log(temp2)/temp2))
              #DF[1,1] <- "1 - (1 - log(changepoint/changepoint)"
              DF[1,1] <- "model 1"
            } #model 1

            else if (k==2) {
              cpiv[temp2,i] <- aveTS*(1 - (temp2-m.bin@cpts[i])/(1 +(temp2-m.bin@cpts[i])))
              DF[2,1] <- "model 2"
            } #model 2

            else if (k==3) {
              cpiv[temp2,i] <- aveTS*(1-(1 /(temp2-m.bin@cpts[i])))
              DF[3,1] <- "model 3"
            } #model 3

            else if (k==4) {
              cpiv[temp2,i] <- aveTS*(1-(1 /(temp2-m.bin@cpts[i])))
              DF[4,1] <- "model 4"
            } #model 4

            else if (k==5) {
              cpiv[temp2,i] <- aveTS*exp(-((temp2-m.bin@cpts[i])^2)/((length(myts)^2)*2.5))
              DF[5,1] <- "model 5"
            } #model 5

            else if (k==6) {
              cpiv[temp2,i] <- aveTS*(1-(1 /(temp2-m.bin@cpts[i])))
              DF[6,1] <- "model 6"
            } #model 6

            else if (k==7) {
              cpiv[temp2,i] <- 1 - (1-(1/log(temp2)))
              DF[7,1] <- "model 7"
            } #model 7

            else if (k==8) {
              cpiv[temp2,i] <- 1
              DF[8,1] <- "model 8"
            } # model 8"

  			    # additional bump functions
            else if (k==9) {
              cpiv[temp2,i] <- StackCutoff (temp2-m.bin@cpts[i], start_level=aveTS, start_transition=1, stop_transition=nSamples-m.bin@cpts[i])
              DF[9,1] <- "model 9"
            } #model 9

            else if (k==10) { # currently use abs to get positive value, can be removed
              cpiv[temp2,i] <- abs(ExpCos (temp2-m.bin@cpts[i], start_point=0, amplitude=2*aveTS, decay.speed=4/nSamples, period=nSamples/3 ))
              DF[10,1] <- "model 10"
            } #model 10

            else if (k==11) { # all standard window functions in signal can be used, here is an example
              # note that the function call is outside the inner loop, here we just take one point at the time
              cpiv[temp2,i] <- aveTS * model11_vector[temp2]
              DF[11,1] <- "model 11"
            } # model 11

            else if (k==12) { # all standard window functions in signal can be used, here is an example
              # note that the function call is outside the inner loop, here we just take one point at the time
              cpiv[temp2,i] <- aveTS * model12_vector[temp2]
              DF[12,1] <- "model 12"
            } # model 12

            else if (k==13) {
              cpiv[temp2,i] <- aveTS * ModulatedExp (temp2-m.bin@cpts[i], start_point=0, amplitude=0.5, decay.speed=4/nSamples, period=nSamples/6 )
              DF[13,1] <- "model 13"
            } #model 13

            else if (k==14) { # sinc based, constrained to be positive, construction insure end point has value 0
              cpiv[temp2,i] <- aveTS * sin (pi*temp2/nSamples) / (pi*temp2/nSamples)
              DF[14,1] <- "model 14"
            } # model 14

            else if (k==15) { # sinc based, constrained to be positive, construction insure end point has value 0 but with more oscillations
              cpiv[temp2,i] <- aveTS * abs(sin (4*pi*temp2/nSamples) / (pi*temp2/nSamples))
              DF[15,1] <- "model 15"
            } # model 15

          }

          if (showThisModel) {
            dev.new() ;
            plot(cpiv[,i], col="blue", type="l", main=sprintf("%s, changepoint %d", DF[k,1], i))
            grid()
          }
        }
        xxx=suppressWarnings(forecast::auto.arima(myts,xreg = cpiv))

        DF[k,3] <- xxx$aic

        DF[[k,4]] <- list(cpiv)
      }

    return(DF)
    }

    suppressWarnings(if (endChangePoint == 0) {
      DF2 <- cpiP2(myts)
    } else {
      for (z in seq(startChangePoint,endChangePoint,step)) {
        if (DF2 == 0) {
          DF2 <- cpiP2(myts,z)
        } else {
          DFtemp <- cpiP2(myts,z)
          DF2 <- rbind(DF2, DFtemp)

        }
      }
    })

    return(DF2)
  }
}
