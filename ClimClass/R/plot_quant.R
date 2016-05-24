NULL
#' @description 'plot' method implementation for 12-month quantile climate charts from output of function \code{\link{thornthwaite}} (Thornthwaite and Mather's water balance).
#' 
#' @param x  a list of quantile data frames of water balance variables to be plotted, as output of function \code{\link{thornthwaite}}.
#' @param save_dir name of destination directory for graphs (if any). Default is \code{NULL}.
#' @param format graphic format of graphs; default is NULL (charts are sent to console).
#' @param variables character vector of variables to be plotted.
#' @param title logic. If \code{TRUE} inserts titles in charts.
#' @param trace_grid logic. If \code{TRUE} (default) adds a grid.
#' @param st_name name to be included into graphs titles. If NULL (default), no title is written.
#' @param l_y_scale_magn magnification of range below lower limit, to set lower y-scale limit; default is 0.1.
#' @param u_y_scale_magn magnification of range above upper limit, to set upper y-scale limit; default is 0.
#' @param leg_pos legend position. Default is "topleft". If NULL, no legend is added.
#' @param ... arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).
#' 
#' @method plot thornthwaite
#' @importFrom graphics legend
#' @export
#' 
#' @title Thornthwaite - Mather's quantile plot
#' 
#' @author Emanuele Eccel
#'
#' @return Charts of quantiles for water balance variables (12-month climatic values). They can be sent to the console or saved as graphic files. 
#'
#' @details
#' Default for plot variables is all those calculated by function thornthwaite: "Precipitation", "Et0", "Storage", "Prec. - Evap.", "Deficit", "Surplus". See function \code{\link{thornthwaite}} for details on variables.
#' 
#' If \code{format} is NULL (default), graphs are sent to the console. Otherwise, a file is produced and saved to the \code{save_dir} directory. Values allowed are: "png", "jpeg", "tiff", "bmp".
#' 
#' \code{l_y_scale_magn} and \code{u_y_scale_magn} are the magnification coefficients (lower and upper, respectively), for y scale. If rng is the range between maximum and minimum values in all sets of series within a plot, the lower limit for y scale will be (rng * \code{l_y_scale_magn}) below the lower value, and the upper limit will be (rng * \code{u_y_scale_magn}) above the upper value of series.
#' 
#' Allowed values for \code{leg_pos} are the same of \code{x} in function \code{\link{legend}}.
#' 
#' Most graphic parameters for functions \code{\link{plot}} and \code{\link{legend}} are accepted.
#' 
#' @note A conflict is generated if parameters already used by the function are passed (e.g. x for \code{\link{legend}}: use \code{leg_pos} instead).
#' 
#' 
#' @examples 
#' 
#' data(Trent_climate)
#' 
#' 
#' # quantiles is the list ("thornthwaite" S3 object)of quantile tables generated 
#' # by function thornthwaite; 
#' # it is the second element of the output list, 
#' # which can be split into two separate lists (see function thornthwaite)
#' sta <- 1     # 1st station in the list of quantile tables
#' q_list=quantiles[[sta]]
#' class(q_list) <- "thornthwaite" ## q_list is coerced to a "thornthwaite" S3 object
#' plot(q_list, 
#' st_name=names(quantiles)[sta], variables=c("Precipitation",  "Et0"), 
#' leg_pos = "topleft", col=c(1:6,1), pch=c(1:6,16),  
#' lty=1, horiz=TRUE,  y.intersp=0.1)
#' 
#' @seealso \code{\link{thornthwaite}}
#'   


plot.thornthwaite<-function(x, save_dir=NULL, format=NULL, variables=c("Precipitation","Et0","Storage","Prec. - Evap.","Deficit","Surplus"), title=TRUE, trace_grid=TRUE, st_name=NULL, u_y_scale_magn=0.2, l_y_scale_magn=0, leg_pos="topleft", ...)
{
  month_names<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  q_list <- x
  for(v in variables)
  {
    d.f.<-as.data.frame(q_list[names(q_list) %in% v])
    y_scale_limits<- c(min(d.f.)-(max(d.f.)-min(d.f.))*l_y_scale_magn ,min(d.f.)+(max(d.f.)-min(d.f.))*(1+u_y_scale_magn))
    legend_text<-substr(row.names(d.f.), 1,(nchar(row.names(d.f.))-1)); legend_text[legend_text== "0"]<-"min";  legend_text[legend_text== "100"]<-"max"
    chart_title<-NULL
    if(!is.null(st_name)) chart_title=paste(st_name, "-", "Monthly", v)
    
    par_l<-list(...)
    if(length(par_l) > 0)
      indices_par_gr<-which(unlist(lapply(X=par_l, FUN= function(x, len){length(x)==len}, len=nrow(d.f.))))
    
    if(!is.null(format))
    {
      if(format=="png") png(filename=paste(save_dir,"/",st_name,"_",v,"_","quant.png", sep="")) else
        if(format=="jpeg") jpeg(filename=paste(save_dir,"/",st_name,"_",v,"_","quant.jpg", sep="")) else
          if(format=="tiff") tiff(filename=paste(save_dir,"/",st_name,"_",v,"_","quant.tif", sep="")) else
            if(format=="bmp") bmp(filename=paste(save_dir,"/",st_name,"_",v,"_","quant.bmp", sep="")) else
              print("Incorrect graphic format", quote=FALSE)
    }
    
    args<-list(x=as.vector(t(d.f.[1,])), main=chart_title, xlab="", ylab="[mm]", ylim=y_scale_limits, type="b", lab=c(12,5,7), xaxt="n")
    par_l_temp<-par_l
    if(length(par_l) > 0)
      par_l_temp[indices_par_gr]<-lapply(X=par_l_temp[indices_par_gr], FUN=function(x, ii){x[ii]}, ii=1)
    
    options(warn=-1)
    
    do.call(what=plot, args=c(args, par_l_temp))
    axis(side=1, labels=month_names, at=1:12)
    for(i in 2:(nrow(d.f.)))
    {
      args<-list(x=as.vector(t(d.f.[i,])), type="b")
      par_l_temp<-par_l
      if(length(par_l) > 0)
        par_l_temp[indices_par_gr]<-lapply(X=par_l_temp[indices_par_gr], FUN=function(x, ii){x[ii]}, ii=i)
      do.call(what=lines, args=c(args, par_l_temp))  
    }
    
    
    if(trace_grid) grid()
    
    if(!is.null(leg_pos))
    {
      par_l_legend<-par_l[names(par_l) %in% names(formals(legend))]
      args_legend<-par_l_legend
      do.call(what=legend, args=c(x=leg_pos, legend= list(legend_text),args_legend ))
    }
    
    if(!is.null(format)) dev.off()
  }
}
