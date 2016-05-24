#' Mass Regression Plot
#'
#' This function takes in a dataframe, the dependent variable, and optionally a character vector of independent variables you want the function to ignore, and produces a regression plot of every variable compared to the dependent variable passed into the function.  It will ignore columns which contain characters and (also optional) factors.
#'
#' @return  Doesn't return a value, per se, but will generate side effects like plotting all the graphs created by the function.  If the save option is used, it will save all generated graphs to a pdf file whose name is specified by the user.
#'
#' @param data data.frame object that contains both the dependent variable and predictor variables you want to plot.
#' 
#' @param dv.var single dependent variable you want to plot your predictors against.
#'
#' @param ignore accepts a character vector of one or more variables you want the function to skip. If nothing is passed through this option, the function will attempt to create a graph plotting the dependent variable and every other column of data.
#'
#' @param save accepts a character.  If the function recieves a character, it will create a pdf file with that name and leave the plots in there.
#' 
#' @param include.factors if TRUE, will also plot factor variables against your dv.  Otherwise it will skip these as regression plots of categorical variables are of imited use.
#' 
#' @param include.se if left TRUE, will shade the area around the regression line with the 95\% confidence interval range.  Setting to FALSE will plot only the regression line to a scatter plot for each paring of variables.
#' 
#' @examples 
#' exam.df <- iris
#' massregplot(exam.df, "Sepal.Length", ignore = "Species")
#' massregplot(exam.df, "Sepal.Length", ignore = c("Species", "Petal.Width"), include.se = FALSE)
#' 
#' @export 

massregplot <- function (data, dv.var, ignore = NULL, save = NULL, include.factors = FALSE, include.se = TRUE) {
  # Things that happen before the function takes place
  mydata <- data
  vars <- colnames(mydata)
  
  # Option to include SE on smooth line or not
  if(include.se == TRUE){
    plot.smooth <- ggplot2::geom_smooth(method = "lm")
  } else { plot.smooth <- ggplot2::geom_smooth(method = "lm", se = FALSE)}
  
  # simple option to save the output. 
  if(class(save) == "NULL"){ 
    } else {
        if(class(save) == "character") {
         grDevices::pdf(save) 
         } else {warning("Save can only accept a character as a filename. Ex: ", save, ".pdf") 
       } 
      }
  
  # Searches through the list of vars and compares them to the list of ignores and the dv.  The function skips those.
  for (i in vars) {
    if (!(i %in% ignore | i == dv.var   )){
      if (include.factors == TRUE) {
        if(!class(mydata[[i]]) %in% "character")
          print(ggplot2::ggplot(mydata, ggplot2::aes(x = get(i), y = get(dv.var))) + ggplot2::geom_point() + plot.smooth + ggplot2::xlab(i) + ggplot2::xlab(dv.var))
      } else if(!class(mydata[[i]]) %in% c("factor", "character")) {
          print(ggplot2::ggplot(mydata, ggplot2::aes(x = get(i), y = get(dv.var))) + ggplot2::geom_point() + plot.smooth + ggplot2::xlab(i) + ggplot2::xlab(dv.var))
      }
    }
  }
  if(class(save) == "character") {
     grDevices::dev.off() 
      for (i in vars) {
        if (!(i %in% ignore | i == dv.var   )){
          if (include.factors == TRUE) {
            if(!class(mydata[[i]]) %in% "character")
              print(ggplot2::ggplot(mydata, ggplot2::aes(x = get(i), y = get(dv.var))) + ggplot2::geom_point() + plot.smooth + ggplot2::xlab(i) + ggplot2::xlab(dv.var))
          } else if(!class(mydata[[i]]) %in% c("factor", "character")) {
            print(ggplot2::ggplot(mydata, ggplot2::aes(x = get(i), y = get(dv.var))) + ggplot2::geom_point() + plot.smooth + ggplot2::xlab(i) + ggplot2::xlab(dv.var))
        }
      }
    }
  }
}
