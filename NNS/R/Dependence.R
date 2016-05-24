#' VN Dependence
#'
#' Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param degree Defaults to 0 for smaller number of observations
#' @param order Number of partial moment quadrants to be generated
#' @param print.map  Displays partition mapping onto plot.  Defaults to TRUE.
#' @keywords dependence
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.dep(x,y)}
#' @export

VN.dep = function( x, y,
                   order=2,
                   degree=0,
                   print.map=TRUE){

  order=order-1
  if(order<1){return("Please Increase the Order Specification")}

  partitioned_df = partition.map(x, y,order)

  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)
  cor.rhos = numeric(0)
  dep.rhos = numeric(0)

  plot(x,y,col='blue',pch=19)



  for(item in unique(partitioned_df$master_part)){


    sub_x = partitioned_df[partitioned_df$master_part == item, 'x']
    sub_y = partitioned_df[partitioned_df$master_part == item, 'y']
    clpm = c(clpm, Co.LPM(degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    cupm = c(cupm, Co.UPM(degree, mean(sub_x),mean(sub_y), sub_x, sub_y))
    dlpm = c(dlpm, D.LPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    dupm = c(dupm, D.UPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))

  if(print.map==TRUE){
    abline(h=mean(y),v=mean(x),lwd=3,col='azure4')
    if(mean(sub_x)<mean(x) && mean(sub_y)<mean(y)){
      segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=1,col='red')
      segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=1,col='red')
      }

    if(mean(sub_x)<mean(x) && mean(sub_y)>mean(y)){
      segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=1,col='red')
      segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=1,col='red')
    }
    if(mean(sub_x)>mean(x) && mean(sub_y)<mean(y)){
      segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=1,col='red')
      segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=1,col='red')
    }
    if(mean(sub_x)>mean(x) && mean(sub_y)>mean(y)){
      segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=1,col='red')
      segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=1,col='red')
    }}

  }


  for(i in 1:(4^order)){

  dep.rhos[i] =  abs((clpm[i]+cupm[i]-dlpm[i]-dupm[i]) / (clpm[i]+cupm[i]+dlpm[i]+dupm[i]))
  }
  for(i in 1:(4^order)){

    cor.rhos[i] =  ((clpm[i]+cupm[i]-dlpm[i]-dupm[i]) / (clpm[i]+cupm[i]+dlpm[i]+dupm[i]))
  }

  cor.rhos = (sum(clpm) +sum(cupm) -sum(dlpm) -sum(dupm))/(sum(clpm)+sum(cupm)+sum(dlpm)+sum(dupm))
  m<- rbind(sum(na.omit(cor.rhos))/length(na.omit(cor.rhos)),
            sum(na.omit(dep.rhos))/length(na.omit(dep.rhos)))

  rownames(m) = c("Correlation","Dependence")
  m[is.nan(m)] <- 0

  #m
  return(m)


}
