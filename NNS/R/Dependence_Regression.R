#' VN Dependence for regression
#'
#' Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.  Routine for regression order, eliminates plotting.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param degree Defaults to 0 for smaller number of observations
#' @param order Number of partial moment quadrants to be generated
#' @keywords dependence
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.dep.reg(x,y)}
#' @export

VN.dep.reg = function( x, y,
                   order=2,
                   degree=0){




  partitioned_df = partition.map(x, y,order)


  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)
  rhos = numeric(0)

  for(item in unique(partitioned_df$master_part)){
    sub_x = partitioned_df[partitioned_df$master_part == item, 'x']
    sub_y = partitioned_df[partitioned_df$master_part == item, 'y']
    clpm = c(clpm, Co.LPM(degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    cupm = c(cupm, Co.UPM(degree, mean(sub_x),mean(sub_y), sub_x, sub_y))
    dlpm = c(dlpm, D.LPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    dupm = c(dupm, D.UPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))



  }


  for(i in 1:order){

    rhos[i] =  abs((clpm[i]+cupm[i]-dlpm[i]-dupm[i]) / (clpm[i]+cupm[i]+dlpm[i]+dupm[i]))
  }



  m<- rbind(VN.cor(x, y,order,degree),sum(na.omit(rhos))/length(na.omit(rhos)))

  rownames(m) = c("Correlation","Dependence")
  m[is.nan(m)] <- 0

  m
  #return(m)

  ### Regression Dependence

  #return(sum(na.omit(rhos))/length(na.omit(rhos)))

}
