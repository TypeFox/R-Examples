NULL
#' @description Plots Bagnouls - Gaussen climatic charts of precipitation and temperature. Conventionally, in this chart the scale of precipitation has a double extension with respect to the scale of temperature (Bagnouls and Gaussen, 1953).
#' 
#' @param clim_norm_sta data frame with climatic normals
#' @param save_dir name of destination directory for graphs (if any).
#' @param format graphical format of graphs; default is NULL.
#' @param main_title main title for all charts; e.g., it may include references to station id. Default is \code{NULL}.
#' @param st_name name to be included into graphs titles. Only for file output. Default is \code{NULL}.
#' @param trace_grid logic. If \code{TRUE} (default) adds a grid.
#' @param tick_step step for Y axis (precipitation). Default is 20 (mm)
#' @param bar_width width of bars in the chart. Default is 30.
#' @param bar_col color of bars. Default is "grey".
#' @param trace_0.line logic. If \code{TRUE} (default), a line at P = 0 and T = 0 is traced.
#' @param ... arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).

#' @title Bagnouls - Gaussen graphs
#' 
#' @author Emanuele Eccel
#'
#' @return Bagnouls - Gaussen's charts of precipitation and temperature. 
#'
#' @details 
#' \code{clim_norm_sta} can be e.g. one element of the output of function \code{\link{climate}}. See \code{examples}.
#' 
#' If \code{format} is NULL (default), graphs are sent to the console. Otherwise, a file is produced and saved. \code{format} is used only if the graphs are to be sent to files. Values allowed are: "png", "jpeg", "tiff", "bmp".
#' 
#' If one or more data are missing, the chart is not processed.
#' 
#' Most graphic parameters for functions \code{\link{plot}}, \code{\link{axis}}, and \code{\link{mtext}} are accepted.
#' 
#' @note A conflict is generated if parameters already used by the function are passed (e.g. \code{col} - use \code{col.main}, \code{col.axis}, ..., instead).
#' @importFrom grDevices png jpeg tiff bmp dev.off
#' @importFrom graphics plot lines axis mtext abline grid
#' @export
#' 
#' @examples
#' 
#' data(Trent_climate)
#' # clima_81_10 can be generated from monthly time series by function \code{\link{climate}}.
#' par(ask=TRUE)
#' for(sta in 1:length(clima_81_10)) {
#'   bagn_gau(clim_norm_sta= clima_81_10 [[sta]], 
#'    main_title=paste(names(clima_81_10[sta]), "  1981-2010")
#' 	, bar_width=40)
#' }
#' 
#' 
#' @seealso \code{\link{climate}}
#'   
#' @references 
#' Bagnouls, F., and Gaussen, H., 1953: Saison seche et indice xerothermique. Docum. pour les Cartes des Prod. Veget. Serie: Generalite, 1 (1953), pp. 1-49


bagn_gau<-function(clim_norm_sta, save_dir=NULL, format=NULL, main_title=NULL, st_name=NULL, trace_grid=TRUE, tick_step=20, bar_width=30, bar_col="grey", trace_0.line=TRUE, ...)
  
{
  
  month_names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  prec<-clim_norm_sta$P
  temp_2<-clim_norm_sta$Tm *2
  
  if(sum(is.na(prec)) > 0 | sum(is.na(temp_2) >0)) print("Missing data not allowed!", quote=FALSE) else
  {
    
    # determine y scale range
    steps_P<-seq(from=-80, to = 1000, by=tick_step)
    
    y_scale_limits_P<-c(0 ,max(prec))
    y_scale_limits_T<-c(min(temp_2) ,max(temp_2))
    y_scale_limits<-c(min(y_scale_limits_P, y_scale_limits_T), max(y_scale_limits_P, y_scale_limits_T))
    y_scale_limits[1]<-max(steps_P[steps_P <= y_scale_limits[1]])
    y_scale_limits[2]<-min(steps_P[steps_P >= y_scale_limits[2]])
    
    if(!is.null(format))
    {
      if(format=="png") png(filename=paste(save_dir,"/",st_name,"_Bagn_Gau.png", sep="")) else
        if(format=="jpeg") jpeg(filename=paste(save_dir,"/",st_name,"_Bagn_Gau.jpg", sep="")) else
          if(format=="tiff") tiff(filename=paste(save_dir,"/",st_name,"_Bagn_Gau.tif", sep="")) else
            if(format=="bmp") bmp(filename=paste(save_dir,"/",st_name,"_Bagn_Gau.bmp", sep="")) else
              print("Incorrect graphic format", quote=FALSE)
    }
    
    
    par_l<-list(...)
    
    args<-list(x=prec, main=main_title, xlab="",ylab="",  ylim=y_scale_limits, type="h", yaxt="n", xaxt="n", lwd=bar_width, lend=1, col=bar_col)
    do.call(what=plot, args=c(args, par_l))
    
    args<-list(x=temp_2, type="l", lwd=2)
    do.call(what=lines, args=c(args, par_l))  
    
    ticks_seq<-seq(from=y_scale_limits[1], to= y_scale_limits[2], by=tick_step)
    
    # x axis 
    args<-list(side=1, at=1:12, labels=month_names)
    do.call(what=axis, args=c(args,par_l))
    # P axis
    at.P<-ticks_seq[ticks_seq >= 0]
    args<-list(side=2, at=at.P, labels=at.P)
    do.call(what=axis, args=c(args,par_l))
    args<-list(side = 2, line = 2, "P [mm]", las=1, adj=0, at=y_scale_limits[2]+ diff(y_scale_limits)*0.15)
    do.call(what=mtext,args=c(args,par_l))
    # T axis
    at.T<-ticks_seq[ticks_seq< max(temp_2)+tick_step]
    args<-list(side=4, at=at.T, labels=at.T/2)
    do.call(what=axis, args=c(args,par_l))
    args<-list(side = 4, line = 2, "T [deg. C]", las=1, adj=1, at=y_scale_limits[2]+ diff(y_scale_limits)*0.15)
    do.call(what=mtext,args=c(args,par_l))
    
    if(trace_0.line) abline(0,0, lty=2) 
    
    if(trace_grid) grid()
    
    if(!is.null(format)) dev.off()
    
  } # else
}
