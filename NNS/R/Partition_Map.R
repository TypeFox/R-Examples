#' Partition Map
#'
#'  Creates partitions based on quadrant means, assigning observations to those quadrants.  Needed for correlation, dependence, regression routines.  Default degree = 1 for area, but routines have their own conditional degree specifications built in
#' @param x Variable 1
#' @param y Variable 2
#' @param order Number of partial moment quadrants to be generated.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{parition.map(x,y)}
#' @export

partition.map = function(x, y,order= ceiling(log10(length(x)))){

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'

  if(order<1){return("Please Increase the Order Specification")}

    order=order-1


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


  return(temp_df[, c('x', 'y', 'master_part')])
}
