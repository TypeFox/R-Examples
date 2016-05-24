#' VN Correlation
#'
#' Returns the nonlinear correlation coefficient based on partial moment quadrants measured by frequency or area.  Degree = 0 is frequency, degree = 1 is area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param degree Defaults to 0 for smaller number of observations
#' @param order Number of partial moment quadrants to be generated
#' @keywords correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.cor(x,y)}
#' @export

VN.cor = function( x, y,order=ceiling(log10(length(x))),
                   degree= ifelse(length(x)<100,0,1)){

  order = order-1

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'


  for(i in 0:(order)){

    for(item in unique(temp_df$master_part)){
      tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
      tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])


      temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 1, sep = '')
      temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 2, sep = '')
      temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 3, sep = '')
      temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 4, sep = '')

    }

    temp_df[,'master_part'] = temp_df[, 'temp_part']
  }

  partitioned_df=temp_df



  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)

  if(order==0){
    clpm = c(clpm, Co.LPM(degree, mean(x),mean(y),x, y))
    cupm = c(cupm, Co.UPM(degree, mean(x),mean(y), x, y))
    dlpm = c(dlpm, D.LPM(degree,degree, mean(x),mean(y),x, y))
    dupm = c(dupm, D.UPM(degree,degree, mean(x),mean(y),x, y))
   }
  else{

    prior.partitioned_df = partitioned_df

    prior.partitioned_df[,'master_part'] = substr(partitioned_df$master_part, 1, nchar(partitioned_df$master_part)-1)

  for(item in unique(prior.partitioned_df$master_part)){

    sub_x = prior.partitioned_df[prior.partitioned_df$master_part == item, 'x']
    sub_y = prior.partitioned_df[prior.partitioned_df$master_part == item, 'y']
    clpm = c(clpm, Co.LPM(degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    cupm = c(cupm, Co.UPM(degree, mean(sub_x),mean(sub_y), sub_x, sub_y))
    dlpm = c(dlpm, D.LPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    dupm = c(dupm, D.UPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
  }}



  nonlin_cor = (sum(clpm) +sum(cupm) -sum(dlpm) -sum(dupm))/(sum(clpm)+sum(cupm)+sum(dlpm)+sum(dupm))

  return(nonlin_cor)


}
