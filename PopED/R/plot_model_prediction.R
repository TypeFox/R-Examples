#' Plot model predictions 
#' 
#' Function plots model predictions for the typical value in the population,
#' individual predictions and data predictions.
#' 
#' @inheritParams RS_opt
#' @inheritParams model_prediction
#' @param separate.groups Should there be separate plots for each group.
#' @param sample.times Should sample times be shown on the plots.
#' @param sample.times.IPRED Should sample times be shown based on the IPRED y-values.
#' @param sample.times.DV Should sample times be shown based on the DV y-values.
#' @param PRED Should a PRED line be drawn.
#' @param IPRED.lines Should IPRED lines be drawn?
#' @param alpha.IPRED.lines What should the transparency for the IPRED.lines be?
#' @param alpha.IPRED What should the tranparency of the IPRED CI?
#' @param sample.times.size What should the size of the sample.times be?
#' @param alpha.DV What should the tranparency of the DV CI?
#' @param DV.lines Should DV lines be drawn?
#' @param DV.points Should DV points be drawn?
#' @param alpha.DV.lines What should the transparency for the DV.lines be?
#' @param alpha.DV.points What should the transparency for the DV.points be?
#' @param sample.times.DV.points TRUE or FALSE.
#' @param sample.times.DV.lines TRUE or FALSE.
#' @param alpha.sample.times.DV.points What should the transparency for the sample.times.DV.points be?
#' @param alpha.sample.times.DV.lines What should the transparency for the sample.times.DV.lines be?
#' @param y_lab The label of the y-axis.
#' @param facet_scales Can be "free", "fixed", "free_x" or "free_y"
#' @param facet_label_names TRUE or FALSE
#' @param IPRED.lines.pctls Should lines be drawn at the chosen percentiles of the IPRED values?  
#' 
#' 
#' @return A \link[ggplot2]{ggplot} object.
#' 
#' @family evaluate_design
#' @family Simulation
#' @family Graphics
#' 
#' @example tests/testthat/examples_fcn_doc/examples_plot_model_prediction.R
#' 
#' @export
#' @import ggplot2
# @import Hmisc
plot_model_prediction <- function(poped.db,
                                  ##models_to_use="all",
                                  model_num_points=100,
                                  ##model_minxt,
                                  ##model_maxxt,
                                  ##groups_to_plot,
                                  ##bShowGraph,bPlotInSame,
                                  separate.groups=F, 
                                  sample.times=T, 
                                  sample.times.IPRED=F,
                                  sample.times.DV=F,
                                  PRED=T,
                                  IPRED=F,
                                  IPRED.lines=F,
                                  IPRED.lines.pctls=F,
                                  alpha.IPRED.lines=0.1,
                                  alpha.IPRED=0.3,
                                  sample.times.size=4,
                                  DV=F,
                                  alpha.DV=0.3,
                                  DV.lines=F,
                                  DV.points=F,
                                  alpha.DV.lines=0.3,
                                  alpha.DV.points=0.3,
                                  sample.times.DV.points=F,
                                  sample.times.DV.lines=F,
                                  alpha.sample.times.DV.points=0.3,
                                  alpha.sample.times.DV.lines=0.3,
                                  y_lab="Model Predictions",
                                  facet_scales="fixed", # could be "free", "fixed", "free_x" or "free_y"
                                  facet_label_names = T, 
                                  ...){
  
  df <-  model_prediction(poped.db,
                          #models_to_use,
                          model_num_points=model_num_points,
                          ##model_minxt,
                          ##model_maxxt,
                          ##groups_to_plot,
                          ...)
  
  if(sample.times){
    df.2 <-  model_prediction(poped.db,
                              ##models_to_use,
                              ##model_num_points=NULL,
                              ##model_minxt,model_maxxt,
                              ##groups_to_plot,
                              ...)
  }
  
  if(IPRED || IPRED.lines || DV || IPRED.lines.pctls){
    df.ipred <-  model_prediction(poped.db,
                                  #models_to_use,
                                  model_num_points=model_num_points,
                                  ##model_minxt,
                                  ##model_maxxt,
                                  ##groups_to_plot,
                                  IPRED=T,
                                  DV=DV,
                                  ...)

    # allow for more experiments to be simulated to get a better idea of the prediction spread
    #     n_experiments=10
    #     df.ipred <- c()
    #     for(i in 1:n_experiments){
    #       df.ipred <-  model_prediction(poped.db,
    #                                     #models_to_use,
    #                                     model_num_points=model_num_points,
    #                                     ##model_minxt,
    #                                     ##model_maxxt,
    #                                     ##groups_to_plot,
    #                                     IPRED=T,
    #                                     DV=DV,
    #                                     ...)
    #     }    
    
    if(sample.times.IPRED || sample.times.DV || sample.times.DV.points || sample.times.DV.lines){
      df.ipred.samples <- df.ipred[df.ipred$Time %in% poped.db$design$xt,]  
    }
  }
  
  ## if((bPlotModel(i))){
  ##     if((bShowGraph)){
  ##         if((!bPlotInSame)){
  ##             fi = figure
  ##         }
  ##         ##hold on
  ##         title(sprintf('Model #d',i))
  
  ##         xlabel('time')
  ##         ylabel('DV')
  ##     }
  
  ## if((bShowGraph)){
  ##     ##hold on
  ##     plot(xt,model_values(i,))
  ##     xlim(matrix(c(minv, maxv),nrow=1,byrow=T))
  ##     #ylim(matrix(c(0,max(model_values(i,))*1.1),nrow=1,byrow=T))
  ##     ##hold off
  ## }
  
  ##df <- df[[1]]
  
  ##library(reshape2)
  ##dat <- melt(df,id=c("Time"))
  ##names(dat) <- c("Time","variable","DV")
  
  many.models <- F
  if(length(unique(df$Model))>1) many.models <- T
  
  many.groups <- F
  if(length(unique(df$Group))>1) many.groups <- T
  
  #library(ggplot2)
  Group <- c()
  Time <- c()
  ID <- c()
   
  if(facet_label_names){
    if(exists("df",inherits=F)) levels(df$Model) <- paste("Model:",levels(df$Model))
    if(exists("df",inherits=F)) levels(df$Group) <- paste("Group:",levels(df$Group))
    if(exists("df.2",inherits=F)) levels(df.2$Model) <- paste("Model:",levels(df.2$Model))
    if(exists("df.2",inherits=F)) levels(df.2$Group) <- paste("Group:",levels(df.2$Group))
    if(exists("df.ipred",inherits=F)) levels(df.ipred$Model) <- paste("Model:",levels(df.ipred$Model))
    if(exists("df.ipred",inherits=F)) levels(df.ipred$Group) <- paste("Group:",levels(df.ipred$Group))
    if(exists("df.ipred.samples",inherits=F)) levels(df.ipred.samples$Model) <- paste("Model:",levels(df.ipred.samples$Model))
    if(exists("df.ipred.samples",inherits=F)) levels(df.ipred.samples$Group) <- paste("Group:",levels(df.ipred.samples$Group))
  }
    
  p <- ggplot(data=df,aes(x=Time,y=PRED)) #+ labs(colour = NULL)
  #p <- ggplot(data=df,aes(x=Time,y=PRED)) #+ labs(colour = NULL)
  if(many.groups) p <- ggplot(data=df,aes(x=Time,y=PRED,group=Group,color=Group,fill=Group)) #+ labs(colour = NULL)
  if(many.models) p <- p + facet_wrap(~Model,scales=facet_scales)
  if(separate.groups && many.models) p <- ggplot(data=df,aes(x=Time,y=PRED,group=Group)) + facet_grid(Model ~ Group,scales=facet_scales)
  if(separate.groups && !many.models) p <- ggplot(data=df,aes(x=Time,y=PRED,group=Group)) + facet_wrap(~Group,scales=facet_scales)
  if(IPRED.lines) p <- p + geom_line(aes(y=IPRED,group=ID),data=df.ipred,alpha=alpha.IPRED.lines)
  if(DV.lines) p <- p + geom_line(aes(y=DV,group=ID),data=df.ipred,alpha=alpha.DV.lines)
  if(DV.points) p <- p + geom_point(aes(y=DV,group=ID),data=df.ipred,alpha=alpha.DV.points)
  if(IPRED) p <- p + stat_summary(data=df.ipred,aes(x=Time,y=IPRED,color=NULL),geom="ribbon",fun.data="median_hilow",alpha=alpha.IPRED)
  if(IPRED.lines.pctls){ 
    p <- p + 
      stat_summary(data=df.ipred,
                   aes(x=Time,y=IPRED),
                   geom="line",
                   fun.y=function(y){quantile(y,probs=0.025)},
                   linetype="dotted")+
      stat_summary(data=df.ipred,
                   aes(x=Time,y=IPRED),
                   geom="line",
                   fun.y=function(y){quantile(y,probs=0.975)},
                   linetype="dotted")
      
  }
  if(DV) p <- p + stat_summary(data=df.ipred,aes(x=Time,y=DV,color=NULL),geom="ribbon",fun.data="median_hilow",alpha=alpha.DV)
  if(PRED) p <- p + geom_line()
  if(sample.times) p <- p+geom_point(data=df.2,size=sample.times.size)#,aes(x=Time,y=DV,group=Group))#,color=Group))
  #if(sample.times) p <- p+geom_point(data=df.2,size=sample.times.size,position = position_jitter(w = 0, h = 10),alpha=0.8)#,aes(x=Time,y=DV,group=Group))#,color=Group))
  #if(sample.times) p <- p+geom_jitter(data=df.2,size=sample.times.size,position = position_jitter(height = 10))#,aes(x=Time,y=DV,group=Group))#,color=Group))
  
  if(sample.times.IPRED) p <- p+stat_summary(data=df.ipred.samples,aes(x=Time,y=IPRED),geom="pointrange",fun.data="median_hilow")
  if(sample.times.DV) p <- p+stat_summary(data=df.ipred.samples,aes(x=Time,y=DV),geom="pointrange",fun.data="median_hilow")
  if(sample.times.DV.lines) p <- p + geom_line(aes(y=DV,group=ID),data=df.ipred.samples,alpha=alpha.sample.times.DV.lines)
  if(sample.times.DV.points) p <- p + geom_point(aes(y=DV,group=ID),data=df.ipred.samples,alpha=alpha.sample.times.DV.points)
  p <- p+ylab(y_lab)
  
  #,scales="free",space="free"
  #ggplot(data=df,aes(x=Time,y=PRED))+  geom_line()+facet_wrap(~Model,scales="free")
  #     p1 + stat_summary(data=df.ipred,aes(x=Time,y=IPRED,fill=Group),geom="ribbon",fun.data="median_hilow",alpha=0.3)
  #     p <- ggplot(data=df.ipred,aes(x=Time,y=IPRED))
  #     p + geom_line(aes(group=ID))
  #     p + geom_line(aes(group=ID,colour=Group))
  #     p + geom_line(alpha=0.3,aes(colour=Group,group=ID))
  #     p + geom_line(aes(group=ID))+facet_grid(~Group) 
  #     
  #     p + geom_line(aes(group=ID)) + stat_summary(geom="ribbon",fun.data="median_hilow",alpha=0.3,fill="red")
  #     
  #     p + stat_summary(geom="ribbon",fun.data="median_hilow",alpha=0.3) +
  #       geom_line(aes(x=Time,y=PRED,group=Group),data=df) + geom_point(aes(x=Time,y=PRED,group=Group),data=df.2,size=4)
  #     
  #     p + stat_summary(geom="ribbon",fun.data="median_hilow",alpha=0.3)+
  #       stat_summary(fun.y = median, #fun.ymin = min, fun.ymax = max,
  #                    colour = "red")
  #     
  #     stat_summary(aes(x=wt,y=obsAUC),geom="ribbon",fun.ymin="min",fun.ymax="max")
  #     p <- ggplot(data=df,aes(x=Time,y=IPRED))
  #     p  + stat_summary(fun.y = median, fun.ymin = min, fun.ymax = max,
  #                       colour = "red")
  #     stat_summary(aes(group=Group,fill=Group),geom="ribbon",fun.data="median_hilow",alpha=0.3) +
  #     fun.ymin="min",fun.ymax="max"
  #     
  #     p+geom_line()+facet_grid(~Group)+geom_point(data=df.2,size=4)
  
  
  return( p ) 
}
