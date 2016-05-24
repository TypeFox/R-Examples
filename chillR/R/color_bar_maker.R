#' Make color scheme for bar plots in outputs of the chillR package
#' 
#' Function to make color schemes for color bar plots in the chillR package.
#' Colors are assigned based on values in two columns of a data frame. One
#' column contains a threshold, below which col3 is assigned. If values are
#' above the threshold, the value in the other column determines the color:
#' col1 if the value is negative, col2 if positive. This function is useful for
#' making the PLS output figures in the chillR package.
#' 
#' 
#' @param column_yn numeric vector containing the data, on which the threshold
#' is to be applied. In the case of the PLS output, this is the data from the
#' VIP column.
#' @param column_quant numeric vector containing the data that determines
#' whether items from column_yn that are above the threshold get assigned col1
#' or col2.
#' @param threshold threshold for values from column_yn to be used for deciding
#' which bars should get col3 and which ones should move on to the next
#' decision step (col1 or col2)
#' @param col1 a color (either a color name or a number) this is applied where
#' column_yn is above the threshold, and column_quant is negative
#' @param col2 a color (either a color name or a number) this is applied where
#' column_yn is above the threshold, and column_quant is positive
#' @param col3 a color (either a color name or a number) this is applied where
#' column_yn is below the threshold
#' @return a vector of colors, which can be used as col argument when making
#' plots
#' @author Eike Luedeling
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' @keywords utility
#' @examples
#' 
#' 
#' PLS_results<-PLS_pheno(
#'   weather_data=make_all_day_table(KA_weather),
#'   split_month=6,   #last month in same year
#'   bio_data=KA_bloom)
#' 
#' colbar<-color_bar_maker(PLS_results$PLS_summary$VIP,PLS_results$PLS_summary$Coef,0.8,
#'        "RED","DARK GREEN","GREY")
#' 
#' @export color_bar_maker
color_bar_maker <-
function(column_yn,column_quant,threshold,col1,col2,col3)
{
  color_bars<-c(rep(NA,length(column_yn)))
  color_bars[which(column_yn>=threshold&column_quant<0)]<-col1
  color_bars[which(column_yn>=threshold&column_quant>=0)]<-col2
  color_bars[which(!column_yn>=threshold)]<-col3
  return(color_bars)}
