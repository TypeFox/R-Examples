#' Observation detection tables
#'
#' Creates a series of tables for dual observer data that shows the number
#' missed and detected for each observer within defined distance classes.
#'
#' @aliases det.tables
#' @param model fitted model from \code{ddf}
#' @param nc number of equal-width bins for histogram
#' @param breaks user define breakpoints
#' @return list object of class "det.tables" \item{Observer1}{table for observer 1}
#' \item{Observer2}{table for observer 2} \item{Duplicates}{histogram counts for duplicates}
#' \item{Pooled}{histogram counts for all observations by either observer}
#' \item{Obs1_2}{table for observer 1 within subset seen by observer 2}
#' \item{Obs2_1}{table for observer 2 within subset seen by observer 1}
#' @author Jeff Laake
#' @export
#' @examples
#' \donttest{
#' data(book.tee.data)
#' region<<-book.tee.data$book.tee.region
#' egdata<<-book.tee.data$book.tee.dataframe
#' samples<<-book.tee.data$book.tee.samples
#' obs<<-book.tee.data$book.tee.obs
#' xx=ddf(mrmodel=~glm(formula=~distance*observer),
#'  dsmodel = ~mcds(key = "hn", formula = ~sex), data = egdata, method = "io",
#' meta.data = list(width = 4))
#' tabs=det.tables(xx,breaks=c(0,.5,1,2,3,4))
#' par(mfrow=c(2,2))
#' plot(tabs,new=FALSE,which=c(1,2,5,6))
#' }
det.tables=function(model,nc=NULL,breaks=NULL){

  # do nothing if we have single observer data
  if(model$method=="ds"){
    stop("This function does not work with single observer data.\n")
  }

  xmat <- process.data(model$data,model$meta.data)$xmat
  xmat$Detected <- factor(xmat$detected,labels=c("Missed","Detected"))
  left <- model$meta.data$left
  width <- model$meta.data$width

  # extract the method
  this.method <- substr(model$method,1,2)

  # Set up default number of classes unless specified
  if(is.null(nc)){
    if(this.method == "io"){
      nc<-round(sqrt(min(
           length(xmat$distance[xmat$observer==1&xmat$detected==1]),
           length(xmat$distance[xmat$observer==2&xmat$detected==1]),
           length(xmat$distance[xmat$observer==1&xmat$timesdetected==2])
          )),0)
    }
    if(this.method == "re"){
      nc <- round(sqrt(
             length(xmat$distance[xmat$observer==2&xmat$detected==1])))
    }
    if(this.method == "tr"){
      nc <- round(sqrt(min(
             length(xmat$distance[xmat$observer==1&xmat$detected==1]),
             length(xmat$distance[xmat$observer==1&xmat$timesdetected==2]) )),0)
    }
  }

  # Set up default break points unless specified
  if(model$meta.data$binned){
    breaks <- model$meta.data$breaks
    nc <- length(breaks)-1
  }else{
    if(is.null(breaks)){
      breaks <- left + ((width-left)/nc)*(0:nc)
    }else{
      nc <- length(breaks)-1
    }
  }

  # Produce tables for each observer, duplicates and pooled
  xmat1 <- xmat[xmat$observer==1,]
  xmat2 <- xmat[xmat$observer==2,]
  obs1 <- with(xmat1,table(cut(distance,breaks,include.lowest=TRUE),Detected))
  obs2 <- with(xmat2,table(cut(distance,breaks,include.lowest=TRUE),Detected))

  if(this.method != "re"){
    obs3 <- with(xmat1[xmat1$timesdetected==2,],
                 table(cut(distance,breaks,include.lowest=TRUE)))
  }else{
    obs3 <- NULL
  }

  if(this.method == "io"){
    obs4 <- with(xmat1,table(cut(distance,breaks,include.lowest=TRUE)))
  }else{
    obs4 <- NULL
  }

  # Produce tables for conditional detection for each observer
  obs1_2 <- as.matrix(with(xmat1[xmat2$detected==1,],
                           table(cut(distance,breaks,include.lowest=TRUE),
                                 Detected)))
  obs1_2 <- cbind(obs1_2,obs1_2[,2]/(obs1_2[,1]+obs1_2[,2]))
  colnames(obs1_2)[3] <- "Prop. detected"
  obs1_2[is.infinite(obs1_2)] <- 0
  if(this.method == "io"){
    obs2_1 <- as.matrix(with(xmat2[xmat1$detected==1,],
                             table(cut(distance,breaks,include.lowest=TRUE),
                                   Detected)))
    obs2_1 <- cbind(obs2_1,obs2_1[,2]/(obs2_1[,1]+obs2_1[,2]))
    colnames(obs2_1)[3] <- "Prop. detected"
    obs2_1[is.infinite(obs2_1)] <- 0
  }else{
    obs2_1 <- NULL
  }

  tab <- list(method=model$method,
              Observer1=obs1,
              Observer2=obs2,
              Duplicates=obs3,
              Pooled=obs4,
              Obs1_2=obs1_2,
              Obs2_1=obs2_1,
              breaks=breaks)

  class(tab) <- "det.tables"
  return(tab)
}
