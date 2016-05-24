###############################################################################
# (c) Copyright IBM Corp. 2015 
# 
# Author: Bradley J. Eck 
###############################################################################


#' Plot Sparkline Table
#' 
#' Generate a table of sparkline plots 
#' 
#' @export 
#' @param df data.frame of values to plot.   
#' @param row.var variable for rows of the table 
#' @param col.vars variables for columns of the table 
#' @param xvar optional name of variable for horizontal axis of sparkline plots
#' @param xrange.labels optional vector of length 2 with labels for the first
#'        and last quantities plotted on x-axis, often a date and/or time
#' 
#' @details Generates a table of 'sparkline' plots of data in df. rows the table correspond to 
#'         different values of row.var. The table's first column gives the value of row.var. The 
#'         remaining columns contain sparkline plots for the values of col.vars.  When xvar is not
#'         provided values are plotted against their index in the extracted vector. The starting
#'         and ending values are labeled. 
#'         Uses layout() function to arrange plots.
#' 
#' @seealso 
#' yaletoolkit and sparkTable packages 
#' @references
#' E. Tufte, Beautiful Evidence, Graphics Press, 2006.
#'  
#' @examples
#' plotSparklineTable( Orange, row.var = 'Tree', col.vars = c('age','circumference'))
#' plotSparklineTable( Loblolly, row.var = 'Seed', col.vars = 'height')
#' ## specify the x variable if you have it, especially if it differs 
#' plotSparklineTable(Theoph, row.var = 'Subject', col.vars = 'conc')
#' ## a warning is normally issued with the ranges of xvar differ 
#' suppressWarnings( plotSparklineTable(Theoph, row.var = 'Subject', col.vars = 'conc', xvar = 'Time'))
plotSparklineTable <- function( df, row.var, col.vars, xvar = NULL, xrange.labels = NULL ){
	
	#this is now a wrapper function 
	slt <- sparklineTable( df, row.var, col.vars, xvar, xrange.labels)
	graphics::plot(slt)
	
} 


