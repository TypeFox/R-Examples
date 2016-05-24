#' Identifiy for a given P&F Table the current vertical price objective 
#' triggered by the last signal reversal.
#' 
#' @param data Input data
#' @param reversal Number of boxes for reversal
#' @param boxsize Size of one box
#' @param log Use logarithmic scale
currentVPOReversalMethod <- function(data,reversal,boxsize,log) {  
  # define local function to identify price objective column
  getPriceObjectiveColumn <- function(data) {
    price.obj.column <- NA
    if (data$column[nrow(data)]-data$column[1]>=3) {
      if (data$status.bs[nrow(data)]=="Buy") {
        column.offset <- 0
        if (data$status.xo[nrow(data)]=="X")
          column.offset <- 1
        columns.to.be.checked <- seq(from=data$column[nrow(data)]-column.offset, to=data$column[1], by=-2)
        for (c in columns.to.be.checked) {
          if (fallingBottom(redData=data,column=c)) {
            price.obj.column <- c+1
            break
          }
        }
      } else if (data$status.bs[nrow(data)]=="Sell") {
        column.offset <- 0
        if (data$status.xo[nrow(data)]=="O")
          column.offset <- 1
        columns.to.be.checked <- seq(from=data$column[nrow(data)]-column.offset, to=data$column[1], by=-2)
        for (c in columns.to.be.checked) {
          if (raisingTop(redData=data,column=c)) {
            price.obj.column <- c+1
            break
          }
        }
      } else {
        stop("Invalid internal status detected!")
      }
    }
    price.obj.column
  }
  
  getPriceObjective <- function(data,min.boxnumber,max.boxnumber,reversal,log) {
    boxnumber <- NA
    price <- NA
    if (data$status.bs[nrow(data)]=="Buy") {
      boxnumber <- min.boxnumber + (max.boxnumber-min.boxnumber+1)*reversal
      # translate price.objective.box into real number
      price <- box2lower(boxnumber=boxnumber,boxsize=boxsize,log=log)
    } else if (data$status.bs[nrow(data)]=="Sell") {
      boxnumber <- max.boxnumber - (max.boxnumber-min.boxnumber+1)*(reversal-1)
      # translate price.objective.box into real number
      price <- box2upper(boxnumber=boxnumber,boxsize=boxsize,log=log)
    } else {
      # should not happen
      stop("Internal Error!")
    }
    list(boxnumber=boxnumber,price=price)
  }
  
  price.objective <- list(boxnumber=NA,price=NA)
  
  ### identify price.obj.column
  price.obj.column <- getPriceObjectiveColumn(data)
  
  if (!is.na(price.obj.column)) {
    # find minimum and maximum boxnumber of price objective column in original data
    min.boxnumber <- NA
    max.boxnumber <- NA
    if (data$status.bs[nrow(data)]=="Buy") {
      max.boxnumber <- max(data$boxnumber[data$column==price.obj.column],na.rm=T)
      min.boxnumber <- min(data$boxnumber[data$column==price.obj.column-1],na.rm=T)+1      
    } else if (data$status.bs[nrow(data)]=="Sell") {
      min.boxnumber <- min(data$boxnumber[data$column==price.obj.column],na.rm=T)
      max.boxnumber <- max(data$boxnumber[data$column==price.obj.column-1],na.rm=T)-1        
    } 
    
    ### determine price objective
    price.objective <- getPriceObjective(data,min.boxnumber,max.boxnumber,reversal,log)
  }
  price.objective
}
