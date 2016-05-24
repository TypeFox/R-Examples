## these methods allow to use the view data button in Rcmdr with reasonable printed output
showData <- function(dataframe,
       colname.bgcolor = "grey50",
       rowname.bgcolor = "grey50",
       body.bgcolor = "white",
       colname.textcolor = "white",
       rowname.textcolor = "white",
       body.textcolor = "black",
       font = "Courier 12",
       maxheight = 30,
       maxwidth = 80,
       title = NULL,
       rowname.bar = "left",
       colname.bar = "top",
       rownumbers = FALSE,
       placement = "-20-40",
       suppress.X11.warnings = TRUE){
  UseMethod("showData")
}

showData.default <- relimp::showData

showData.design <- function(dataframe, colname.bgcolor = "grey50",
       rowname.bgcolor = "grey50",
       body.bgcolor = "white",
       colname.textcolor = "white",
       rowname.textcolor = "white",
       body.textcolor = "black",
       font = "Courier 12",
       maxheight = 30,
       maxwidth = 80,
       title = NULL,
       rowname.bar = "left",
       colname.bar = "top",
       rownumbers = FALSE,
       placement = "-20-40",
       suppress.X11.warnings = TRUE) {
   datnam <- deparse(substitute(dataframe))
   if (!"design" %in% class(dataframe))
       stop("This method is for class design data frames only.")
   showData(undesign(dataframe),colname.bgcolor=colname.bgcolor,
      rowname.bgcolor=rowname.bgcolor, body.bgcolor=body.bgcolor,
      colname.textcolor=colname.textcolor, rowname.textcolor=rowname.textcolor,
      body.textcolor=body.textcolor, font=font,maxheight=maxheight,
      maxwidth=maxwidth,title=datnam, rowname.bar=rowname.bar, colname.bar=colname.bar,
      rownumbers=rownumbers,placement=placement,suppress.X11.warnings=suppress.X11.warnings)
   }
