#' FSD
#'
#' Uni-directional test of first degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x variable
#' @param y variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{FSD(x,y)}
#' @export

FSD <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = numeric(0)
  LPM_y_sort = numeric(0)

  output_x <- vector("numeric", length(x))
  output_y <- vector("numeric", length(x))

  if(min(y)>=min(x)) {return(0)} else {


  for (i in 1:length(Combined)){

    ## Indicator function ***for all values of x and y*** as the CDF target
    if(LPM(0,Combined_sort[i],y)-LPM(0,Combined_sort[i],x)>=0 )
    {output_x[i]<-0} else { break } }


  for (i in 1:length(Combined)){
    if(LPM(0,Combined_sort[i],x)-LPM(0,Combined_sort[i],y)>=0 )
    {output_y[i]<-0} else { break }

  }


  for (j in 1:length(Combined_sort)){
    LPM_x_sort[j] = LPM(0,Combined_sort[j],x)
    LPM_y_sort[j] = LPM(0,Combined_sort[j],y)
  }

  ifelse(length(output_x)==length(Combined) & min(x)>=min(y),return(1),return(0))

  }
}

#' VN.SSD.uni
#'
#' Uni-directional test of second degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x variable
#' @param y variable
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.SSD.uni(x,y)}
#' @export


VN.SSD.uni <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = numeric(0)
  LPM_y_sort = numeric(0)

  output_x <- vector("numeric", length(x))
  output_y <- vector("numeric", length(x))

  if(min(y)>=min(x)) {return(0)} else {


    for (i in 1:length(Combined)){

      ## Indicator function ***for all values of x and y*** as the CDF target
      if(LPM(1,Combined_sort[i],y)-LPM(1,Combined_sort[i],x)>=0 )
      {output_x[i]<-0} else { break } }


    for (i in 1:length(Combined)){
      if(LPM(1,Combined_sort[i],x)-LPM(1,Combined_sort[i],y)>=0 )
      {output_y[i]<-0} else { break }

    }


    for (j in 1:length(Combined_sort)){
      LPM_x_sort[j] = LPM(1,Combined_sort[j],x)
      LPM_y_sort[j] = LPM(1,Combined_sort[j],y)
    }

    ifelse(length(output_x)==length(Combined) & min(x)>=min(y),return(1),return(0))

  }
}


#' TSD
#'
#' Uni-directional test of third degree stochastic dominance using lower partial moments used in SD Efficient Set routine.
#' @param x variable
#' @param y variable
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{TSD(x,y)}
#' @export


TSD <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = numeric(0)
  LPM_y_sort = numeric(0)

  output_x <- vector("numeric", length(x))
  output_y <- vector("numeric", length(x))

  if(min(y)>=min(x) | mean(y)>=mean(x)) {return(0)} else {


    for (i in 1:length(Combined)){

      ## Indicator function ***for all values of x and y*** as the CDF target
      if(LPM(2,Combined_sort[i],y)-LPM(2,Combined_sort[i],x)>=0 )
      {output_x[i]<-0} else { break } }


    for (i in 1:length(Combined)){
      if(LPM(2,Combined_sort[i],x)-LPM(2,Combined_sort[i],y)>=0 )
      {output_y[i]<-0} else { break }

    }


    for (j in 1:length(Combined_sort)){
      LPM_x_sort[j] = LPM(2,Combined_sort[j],x)
      LPM_y_sort[j] = LPM(2,Combined_sort[j],y)
    }

    ifelse(length(output_x)==length(Combined) & min(x)>=min(y),return(1),return(0))

  }
}
