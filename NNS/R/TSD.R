#' VN TSD Test
#'
#' Bi-directional test of third degree stochastic dominance using lower partial moments.
#' @param x variable
#' @param y variable
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2016) "LPM Density Functions for the Computation of the SD Efficient Set." Journal of Mathematical Finance, 6, 105-126. \url{http://www.scirp.org/Journal/PaperInformation.aspx?PaperID=63817}.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.TSD(x,y)}
#' @export

VN.TSD <- function(x,y){

  x_sort <- sort(x, decreasing=FALSE)
  y_sort <- sort(y, decreasing=FALSE)

  Combined = c(x_sort,y_sort)
  Combined_sort = sort(Combined, decreasing=FALSE)

  LPM_x_sort = numeric(0)
  LPM_y_sort = numeric(0)

  output_x <- vector("numeric", length(x))
  output_y <- vector("numeric", length(x))


  for (i in 1:length(Combined)){

    if(LPM(2,Combined_sort[i],y)-LPM(2,Combined_sort[i],x)>=0 )
    {output_x[i]<-0} else { break }}

  for (i in 1:length(Combined)){

    if(LPM(2,Combined_sort[i],x)-LPM(2,Combined_sort[i],y)>=0 )
    {output_y[i]<-0} else { break }

  }



  for (j in 1:length(Combined_sort)){
    LPM_x_sort[j] = LPM(2,Combined_sort[j],x)
    LPM_y_sort[j] = LPM(2,Combined_sort[j],y)
  }

  plot(LPM_x_sort, type = "l", lwd =3,col = "red", main = "TSD", ylab = "Area of Cumulative Distribution",
       ylim = c(min(c(LPM_y_sort,LPM_x_sort)),max(c(LPM_y_sort,LPM_x_sort))))
  lines(LPM_y_sort, type = "l", lwd =3,col = "blue")
  legend("topleft", c("X","Y"), lwd=10,
         col=c("red","blue"))

   ifelse (length(output_x)==length(Combined) & min(x)>=min(y) & mean(x)>=mean(y),"X TSD Y",
          ifelse (length(output_y)==length(Combined) & min(y)>=min(x)& mean(y)>=mean(x),"Y TSD X","NO TSD EXISTS"))
}

