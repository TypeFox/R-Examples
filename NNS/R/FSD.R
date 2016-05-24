#' VN FSD Test
#'
#' Bi-directional test of first degree stochastic dominance using lower partial moments.
#' @param x variable
#' @param y variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.FSD(x,y)}
#' @export



VN.FSD <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = numeric(0)
  LPM_y_sort = numeric(0)

  output_x <- vector("numeric", length(x))
  output_y <- vector("numeric", length(x))


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

    plot(LPM_x_sort, type = "l", lwd =3,col = "red", main = "FSD", ylab = "Probability of Cumulative Distribution")
    lines(LPM_y_sort, type = "l", lwd =3,col = "blue")
    legend("topleft", c("X","Y"), lwd=10,
           col=c("red","blue"))

     ## Verification of ***0 instances*** of CDFx > CDFy, and conversely of CDFy > CDFx
    ifelse (length(output_x)==length(Combined) & min(x)>=min(y),"X FSD Y",
           ifelse (length(output_y)==length(Combined) & min(y)>=min(x),"Y FSD X","NO FSD EXISTS"))
}

