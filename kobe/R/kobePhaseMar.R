utils::globalVariables(c("stock","..count..","group","harvest","run"))
#' @importFrom grid grid.newpage pushViewport viewport grid.layout unit
#' 
#' @name kobePhaseMar
#' @title Kobe Phase Plot with Marginal Densities
#' @description 
#' Plots a Kobe Phase Plot with marginal densities
#' '      
#' @aliases kobePhaseMar
#' 
#' @param  pts data.frame with stock and harvest variables respresenting points for a year
#' @param  trks data.frame with stock and harvest variables respresenting tracks, NULL by default
#' @param  mns logical, do you want mean by year of pts be plotted, FALSE  by default
#' @param  size of points1  by default
#' @param  xlab character or expression for x axis title, expression(B:B[MSY])  by default
#' @param  ylab character or expression for x axis title, expression(F:F[MSY])  by default
#' @param  xlim upper limit of x axis, 2  by default
#' @param  ylim upper limit of x axis, 2  by default
#' @param  col colorRampPalette(c("orange","blue"),space="Lab")  by default
#' @param  shade .5  by default
#' @param  col2 grey(shade)  by default
#' @param  col3 grey(shade*1.1)  by default
#' 
#' @return a list with 3 ggplot objects, the 2 marginal densities and the phase plot
#' @export kobePhaseMar
#' @docType methods
#' @rdname  kobePhaseMar-method
#' 
#' @examples
#'    \dontrun{rnorm(10)}
kobePhaseMar=function(pts,trks=NULL,mns=FALSE,size=1,
             xlab=expression(B:B[MSY]),
             ylab=expression(F:F[MSY]),
             xlim=2,ylim=xlim,
             col =colorRampPalette(c("orange","blue"),space="Lab"),
             shade=.5,col2=grey(shade),col3=grey(shade*1.1)){
     
    if (!("run" %in% names(pts)))
       pts=cbind(pts,run=factor(1))
    if (!is.null(trks) & !("run" %in% names(trks)))
      trks=cbind(trks,run=factor(1))
 
    if ("function" %in% is(col))
       col=col(length(unique(pts$run)))
   
    if (length(size)==1) size=rep(size,2)
    
    ##### Density plots   #############################################################################################
    # stock density plot
    dS<-ggplot(pts) + 
          geom_density(aes(x = stock, y =  ..count.., group=run), fill=col2, col=col3, position = "stack") + 
          geom_density(aes(x = stock, y = -..count.., fill =run, alpha=0.4)) + 
          geom_vline(xintercept=1,col="red")       +
              scale_x_continuous(limits=c(0,xlim)) +
              scale_fill_manual(values=col)        +
              xlab(xlab) + ylab(ylab)              +
              theme(legend.position = "none", 
                    axis.title.y = element_text(colour='NA'), 
                    axis.text.y  = element_text(colour="NA", angle=90), 
                    axis.ticks.y = element_line(colour="NA"),
                    axis.ticks =   element_line(colour="NA"),
                    
                    axis.title.x = element_blank(), 
                    axis.text.x  = element_blank(), 
                    axis.ticks.x = element_blank(), 
                    
                    plot.margin = unit(c(0, 0, 0, 1), "lines"),
                    panel.background = element_rect(fill   ="NA", colour="NA"), 
                    panel.border     = element_rect(fill   ="NA", colour="NA"), 
                    panel.grid.major = element_line(colour ="NA"), 
                    panel.grid.minor = element_line(colour ="NA")              )
    
    # second density plot, oriented vertically (hence the 'coord_flip()' at the end
      dH<-ggplot(pts) + 
            geom_density(aes(x = harvest, y =  ..count.., group=run), fill=col2, col=col3, position = "stack") + 
            geom_density(aes(x = harvest, y = -..count..,               fill=run, alpha=0.4)) + 
            geom_vline(xintercept=1,col="red")  +
                scale_x_continuous(limits=c(0,ylim))   +
                scale_fill_manual(values=col)          +
                xlab(xlab) + ylab(ylab)                +
                theme(legend.position = "none", 
                      axis.title.x = element_text(colour ='NA'), 
                      axis.text.x  = element_text(colour ="NA"), 
                      axis.ticks.x = element_line(colour ="NA"),
                      axis.ticks =   element_line(colour ="NA"),
                      
                      axis.title.y = element_blank(), 
                      axis.text.y  = element_blank(), 
                      axis.ticks.y = element_blank(), 
                      
                      plot.margin = unit(c(0, 0, 1, 0), "lines"),
                      panel.background = element_rect(fill   ="NA", colour ="NA"), 
                      panel.border     = element_rect(fill   ="NA", colour ="NA"), 
                      panel.grid.major = element_line(colour ="NA"), 
                      panel.grid.minor = element_line(colour ="NA")                    
                )
  
    # kobe phase plot
    kC=kobePhase(pts) +
       geom_point(aes(stock,harvest,col=run,group=run),size=size[1]) +  
       scale_y_continuous(limits=c(0,ylim)) +
       scale_x_continuous(limits=c(0,xlim)) +
       scale_colour_manual(values=col)      +
       xlab(xlab) + ylab(ylab)              +
       theme(legend.position = "none",
             axis.text.y  = element_text(colour="grey", angle=90), 
             plot.margin = unit(c(0, 0, 1, 1), "lines")
       )
   
#     if (length(run)>1){
#         dS=dS+scale_fill_manual(values=col)
#         dH=dH+scale_fill_manual(values=col)
#         kC=kC+scale_colour_manual(values=col)      
#         }
#     
    if (mns)
       kC=kC+geom_point(aes(stock,harvest,col=run,group=run),size=6.0*size[1], colour="black",  data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest)))) +
             geom_point(aes(stock,harvest,col=run,group=run),size=4.5*size[1], colour="cyan",   data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest))))
    if (!is.null(trks))
       kC=kC+geom_path(aes(stock,harvest, col=run,group=run),size=1*size[2], data=trks)   
    
#     fnVP=function(dH,dS,kC){
#         vplayout <- function(x, y)
#         viewport(layout.pos.row = x, layout.pos.col = y)
#           
#         grid.newpage()
#         pushViewport(viewport(layout = grid.layout(5, 5)))  # 5 by 5 grid
#         print(dS, vp=vplayout(1,1:4))                       # the first density plot will occupy the top of the grid
#         print(dH +coord_flip(), vp=vplayout(2:5,5))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
#         print(kC, vp=vplayout(2:5,1:4))                     # the main x/y plot will instead spread across most of the grid
#         }
#     
#     fnVP(dH,dS,kC)
    print(kC)

    invisible(list(harvest=NULL,stock=NULL,phase=kC))}


kobePhaseMar2=function(pts,trks=NULL,mns=FALSE,size=1,
                       xlab=expression(B:B[MSY]),
                       ylab=expression(F:F[MSY]),
                       xlim=2,ylim=xlim,
                       col =colorRampPalette(c("orange","blue"),space="Lab"),
                       shade=.5,col2=grey(shade),col3=grey(shade*1.1)){
  
  if (!("run" %in% names(pts)))
    pts=cbind(pts,run=factor(1))
  if (!is.null(trks) & !("run" %in% names(trks)))
    trks=cbind(trks,run=factor(1))
  
  if ("function" %in% is(col))
    col=col(length(unique(pts$run)))
  
  if (length(size)==1) size=rep(size,2)
  
  ##### Density plots   #############################################################################################
  # stock density plot
  dS<-ggplot(pts) + 
    #geom_density(aes(x = stock, y =  ..count.., group=run), fill=col2, col=col3, position = "stack") + 
    geom_density(aes(x = stock, y = ..count.., fill =run, alpha=0.4)) + 
    geom_vline(xintercept=1,col="red")       +
    scale_x_continuous(limits=c(0,xlim)) +
    scale_fill_manual(values=col)        +
    xlab(xlab) + ylab(ylab)              +
    theme(legend.position = "none", 
          axis.title.y = element_text(colour='NA'), 
          axis.text.y  = element_text(colour="NA", angle=90), 
          axis.ticks.y = element_line(colour="NA"),
          axis.ticks =   element_line(colour="NA"),
          
          axis.title.x = element_blank(), 
          axis.text.x  = element_blank(), 
          axis.ticks.x = element_blank(), 
          
          plot.margin = unit(c(0, 0, 0, 1), "lines"),
          panel.background = element_rect(fill   ="NA", colour="NA"), 
          panel.border     = element_rect(fill   ="NA", colour="NA"), 
          panel.grid.major = element_line(colour ="NA"), 
          panel.grid.minor = element_line(colour ="NA")              )
  
  # second density plot, oriented vertically (hence the 'coord_flip()' at the end
  dH<-ggplot(pts) + 
    #geom_density(aes(x = harvest, y =  ..count.., group=run), fill=col2, col=col3, position = "stack") + 
    geom_density(aes(x = harvest, y = ..count..,               fill=run, alpha=0.4)) + 
    geom_vline(xintercept=1,col="red")  +
    scale_x_continuous(limits=c(0,ylim))   +
    scale_fill_manual(values=col)          +
    xlab(xlab) + ylab(ylab)                +
    theme(legend.position = "none", 
          axis.title.x = element_text(colour ='NA'), 
          axis.text.x  = element_text(colour ="NA"), 
          axis.ticks.x = element_line(colour ="NA"),
          axis.ticks =   element_line(colour ="NA"),
          
          axis.title.y = element_blank(), 
          axis.text.y  = element_blank(), 
          axis.ticks.y = element_blank(), 
          
          plot.margin = unit(c(0, 0, 1, 0), "lines"),
          panel.background = element_rect(fill   ="NA", colour ="NA"), 
          panel.border     = element_rect(fill   ="NA", colour ="NA"), 
          panel.grid.major = element_line(colour ="NA"), 
          panel.grid.minor = element_line(colour ="NA")                    
    )
  
  # kobe phase plot
  kC=kobePhase(pts) +
    geom_point(aes(stock,harvest,col=run,group=run),size=size[1]) +  
    scale_y_continuous(limits=c(0,ylim)) +
    scale_x_continuous(limits=c(0,xlim)) +
    scale_colour_manual(values=col)      +
    xlab(xlab) + ylab(ylab)              +
    theme(legend.position = "none",
          axis.text.y  = element_text(colour="grey", angle=90), 
          plot.margin = unit(c(0, 0, 1, 1), "lines")
    )
  
  #     if (length(run)>1){
  #         dS=dS+scale_fill_manual(values=col)
  #         dH=dH+scale_fill_manual(values=col)
  #         kC=kC+scale_colour_manual(values=col)      
  #         }
  #     
  if (mns)
    kC=kC+geom_point(aes(stock,harvest,col=run,group=run),size=6.0*size[1], colour="black",  data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest)))) +
    geom_point(aes(stock,harvest,col=run,group=run),size=4.5*size[1], colour="cyan",   data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest))))
  if (!is.null(trks))
    kC=kC+geom_path(aes(stock,harvest, col=run,group=run),size=1*size[2], data=trks)   
  
  fnVP=function(dH,dS,kC){
    vplayout <- function(x, y)
      viewport(layout.pos.row = x, layout.pos.col = y)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5, 5)))  # 5 by 5 grid
    print(dS, vp=vplayout(1,1:4))                       # the first density plot will occupy the top of the grid
    print(dH +coord_flip(), vp=vplayout(2:5,5))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
    print(kC, vp=vplayout(2:5,1:4))                     # the main x/y plot will instead spread across most of the grid
  }
  
  fnVP(dH,dS,kC)
  
  invisible(list(harvest=dH,stock=dS,phase=kC))}

kobePhaseMar3=function(pts,trks=NULL,mns=FALSE,size=1,
                       xlab=expression(B:B[MSY]),
                       ylab=expression(F:F[MSY]),
                       xlim=2,ylim=xlim,
                       col =colorRampPalette(c("orange","blue"),space="Lab"),
                       shade=.5,col2=grey(shade),col3=grey(shade*1.1)){
  
  if (!("run" %in% names(pts)))
    pts=cbind(pts,run=factor(1))
  if (!is.null(trks) & !("run" %in% names(trks)))
    trks=cbind(trks,run=factor(1))
  
  if ("function" %in% is(col))
    col=col(length(unique(pts$run)))
  
  if (length(size)==1) size=rep(size,2)

  ##### Density plots   #############################################################################################
  # stock density plot
  dS<-ggplot(pts) + 
    geom_density(aes(x = stock, y =  ..count.., group=run), fill=col2, col=col3, position = "stack") + 
    #geom_density(aes(x = stock, y = ..count.., fill =run, alpha=0.4)) + 
    geom_vline(xintercept=1,col="red")       +
    scale_x_continuous(limits=c(0,xlim)) +
    scale_fill_manual(values=col)        +
    xlab(xlab) + ylab(ylab)              +
    theme(legend.position = "none", 
          axis.title.y = element_text(colour='NA'), 
          axis.text.y  = element_text(colour="NA", angle=90), 
          axis.ticks.y = element_line(colour="NA"),
          axis.ticks =   element_line(colour="NA"),
          
          axis.title.x = element_blank(), 
          axis.text.x  = element_blank(), 
          axis.ticks.x = element_blank(), 
          
          plot.margin = unit(c(0, 0, 0, 1), "lines"),
          panel.background = element_rect(fill   ="NA", colour="NA"), 
          panel.border     = element_rect(fill   ="NA", colour="NA"), 
          panel.grid.major = element_line(colour ="NA"), 
          panel.grid.minor = element_line(colour ="NA")              )
  
  # second density plot, oriented vertically (hence the 'coord_flip()' at the end
  dH<-ggplot(pts) + 
    geom_density(aes(x = harvest, y =  ..count.., group=run), fill=col2, col=col3, position = "stack") + 
    #geom_density(aes(x = harvest, y = ..count..,               fill=run, alpha=0.4)) + 
    geom_vline(xintercept=1,col="red")  +
    scale_x_continuous(limits=c(0,ylim))   +
    scale_fill_manual(values=col)          +
    xlab(xlab) + ylab(ylab)                +
    theme(legend.position = "none", 
          axis.title.x = element_text(colour ='NA'), 
          axis.text.x  = element_text(colour ="NA"), 
          axis.ticks.x = element_line(colour ="NA"),
          axis.ticks =   element_line(colour ="NA"),
          
          axis.title.y = element_blank(), 
          axis.text.y  = element_blank(), 
          axis.ticks.y = element_blank(), 
          
          plot.margin = unit(c(0, 0, 1, 0), "lines"),
          panel.background = element_rect(fill   ="NA", colour ="NA"), 
          panel.border     = element_rect(fill   ="NA", colour ="NA"), 
          panel.grid.major = element_line(colour ="NA"), 
          panel.grid.minor = element_line(colour ="NA")                    
    )
  
  # kobe phase plot
  kC=kobePhase(pts) +
    geom_point(aes(stock,harvest,col=run,group=run),size=size[1]) +  
    scale_y_continuous(limits=c(0,ylim)) +
    scale_x_continuous(limits=c(0,xlim)) +
    scale_colour_manual(values=col)      +
    xlab(xlab) + ylab(ylab)              +
    theme(legend.position = "none",
          axis.text.y  = element_text(colour="grey", angle=90), 
          plot.margin = unit(c(0, 0, 1, 1), "lines")
    )
  
  #     if (length(run)>1){
  #         dS=dS+scale_fill_manual(values=col)
  #         dH=dH+scale_fill_manual(values=col)
  #         kC=kC+scale_colour_manual(values=col)      
  #         }
  #     
  if (mns)
    kC=kC+geom_point(aes(stock,harvest,col=run,group=run),size=6.0*size[1], colour="black",  data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest)))) +
    geom_point(aes(stock,harvest,col=run,group=run),size=4.5*size[1], colour="cyan",   data=ddply(pts,.(run),function(x) data.frame(stock=median(x$stock),harvest=median(x$harvest))))
  if (!is.null(trks))
    kC=kC+geom_path(aes(stock,harvest,col=run,group=run),size=1*size[2], data=trks)   
  
  fnVP=function(dH,dS,kC){
    vplayout <- function(x, y)
      viewport(layout.pos.row = x, layout.pos.col = y)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5, 5)))  # 5 by 5 grid
    print(dS, vp=vplayout(1,1:4))                       # the first density plot will occupy the top of the grid
    print(dH +coord_flip(), vp=vplayout(2:5,5))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
    print(kC, vp=vplayout(2:5,1:4))                     # the main x/y plot will instead spread across most of the grid
  }
  
  fnVP(dH,dS,kC)
  
  invisible(list(harvest=dH,stock=dS,phase=kC))}

