#in development code
#[TBC - NUMBER] functions 

#panel.polarPlot
#panel.polarFrame
#panel.polarAxes
#panel.polarGrid
#panel.polarLabels

#NOTE: much borrowed from... 

#to do

############################
#repairs
############################

#grid and labels need a fix re rlim
##should it go (0, max)?

#axes not showing?



##############################
##############################
##panel.polarPlot
##############################
##############################

panel.polarPlot <- function(x = NULL, y = NULL, r = NULL, theta = NULL, ..., 
         data.panel = panel.loaPlot,  
         loa.settings = FALSE, plot = TRUE, process = TRUE){

#################
#panel to link 
#cex and col to 
#colorkey and z
#################

#think about
##################
#update from list(output="col") for colHandler?
#to make more robust

#currently cannot pass any formals
#for loaplot 
##think about this



    ###################
    #return safe.mode info
    ###################

##############
#aspect needs thinking about
##############

    if(loa.settings)
        return(list(group.args= c("col"),
                    zcase.args= c("pch"),
                    default.settings = list(local.scales = TRUE,
                                            local.scales.panel = panel.polarFrame,
                                            grid = TRUE, axes = TRUE, labels = TRUE,
                                            allowed.scales =c("r", "theta"),  
                                            disallowed.scales = c("x", "y"),
                                            aspect = "loa.iso", 
                                            reset.xylims = c("refit.xylims", "max.xylims"),
                                            load.lists = c("grid", "axes", "labels"),                                            
                                            key.fun = "draw.loaPlotZKey")))

    if(process){

#####################
#should this be if r and theta NULL
#####################
        
        r <- x
        theta <- y
        x <- r * sin(pi * theta/180)
        y <- r * cos(pi * theta/180) 
        if(!plot) return(list(x = x, y = y, r = r, theta = theta))
    }

    if(plot){

        extra.args <- listUpdate(list(...), 
                                 list(x = x, y = y, r = r, theta = theta, 
                                      plot = plot, process = process))

###################
#might want to allow user to set this
#might want to be able to turn this off
#so not local.scales.panels
###################

        if(extra.args$local.scales)
            if(is.function(extra.args$local.scales.panel)){
                do.call(extra.args$local.scales.panel, extra.args)
                extra.args$grid <- NULL
            }
        do.call(data.panel, extra.args)
    }
}


###########################
###########################
##panel.polarFRame
##panel.panelAxes
##panel.polarGrid
##panel.polarLabels
###########################
###########################


panel.polarFrame <- function(..., grid = TRUE, axes = TRUE, labels = TRUE, 
         panel.scales = NULL,
         grid.panel = panel.polarGrid, 
         axes.panel = panel.polarAxes, 
         labels.panel = panel.polarLabels
){

    extra.args <- list(...)

####################################
#might need to rethink this
#temp fix so grid has range like 
#default
####################################

    if(!is.null(extra.args$rlim))
        extra.args$rlim <- c(0, extra.args$rlim)

#####################################
#might need to make these next calls
#update from extra args
#####################################


    if(isGood4LOA(grid))
           grid.panel(rlim=extra.args$rlim, grid.theta=(0:12)*30, 
                       panel.scales = panel.scales, grid = grid) 
           
    if(isGood4LOA(axes))
           axes.panel(rlim=extra.args$rlim, axes.theta=(0:3)*90, 
                      panel.scales = panel.scales, axes = axes) 

    if(isGood4LOA(labels))
           labels.panel(rlim=extra.args$rlim, labels.theta=((0:3)*90)+45, 
                        panel.scales = panel.scales, labels = labels) 

}


##############################################
##############################################

panel.polarAxes <- function(axes.theta = NULL, axes.r = NULL,
         thetalim = NULL, rlim = NULL, ..., 
         axes = NULL, panel.scales = NULL){

    extra.args <- list(...)

    if (!is.list(panel.scales)) panel.scales <- list()
    if (!is.list(axes)) axes <- list()

    panel.scales <- listUpdate(list(draw = TRUE, arrows = FALSE, tick.number = 5, 
                                    abbreviate = FALSE, minlength = 4, tck = 1, 
                                    col.line = 1, cex = 0.8), 
                              panel.scales)


    theta.par <- getPlotArgs("axis.line", local.resets = panel.scales,
                             user.resets = axes, elements = "theta",
                             defaults.only = FALSE)


    if(!is.null(theta.par$theta))
        axes.theta <- theta.par$theta
    if(!is.null(theta.par$grid.theta))
        axes.theta <- theta.par$grid.theta
    if(is.null(axes.theta)){
        axes.theta <- if(!is.null(theta.par$at))
                           theta.par$at else
                               if(!is.null(thetalim))
                                   pretty(thetalim, theta.par$tick.number) else
                                      (0:12) * 30    
    }


    r.par <- getPlotArgs("axis.line", local.resets = panel.scales,
                         user.resets = axes, elements = "r", 
                         defaults.only=FALSE)

    if(!is.null(r.par$r))
        axes.r <- r.par$r
    if(!is.null(r.par$grid.r))
        axes.r <- r.par$grid.r
    if(is.null(axes.r)){
        axes.r <- if(!is.null(r.par$at))
                        r.par$at else
                            if(!is.null(rlim))
                                pretty(rlim, r.par$tick.number) else
                                      NULL    
    }



#    for(i in grid.r){
#        theta <- seq(min(grid.theta), max(grid.theta), 
#                     length.out = 10 * (max(grid.theta)-min(grid.theta)))
#        r <- rep(i, length(theta))
#        x <- r * sin(pi * theta/180)
#        y <- r * cos(pi * theta/180)
#        temp <- listUpdate(list(x=x, y=y, col="green"), theta.par)
#        temp$col.line <- temp$col
#        do.call(panel.lines, temp)
#    }

##mix fix
##same in polar.polarGrid


    for(i in axes.theta){
        theta <- c(0, i)
#        r <- c(min(axes.r), max(axes.r))
        r <- c(0, max(axes.r))
        x <- r * sin(pi * theta/180)
        y <- r * cos(pi * theta/180)
        temp <- listUpdate(list(x=x, y=y), r.par)
        temp$col.line <- temp$col
        do.call(panel.lines, temp)
    }

}









##############################################
##############################################


panel.polarGrid <- function(grid.theta = NULL, grid.r = NULL,
         thetalim = NULL, rlim = NULL, ..., 
         grid = NULL, panel.scales = NULL){

######################
#might want to rethink this
#for grid being the ....
#then grid.theta could be dropped


    extra.args <- list(...)

    if (!is.list(panel.scales)) panel.scales <- list()
    if (!is.list(grid)) grid <- list()

    panel.scales <- listUpdate(list(draw = TRUE, arrows = FALSE, tick.number = 5, 
                                    abbreviate = FALSE, minlength = 4, tck = 1, 
                                    col = "lightgrey", col.line = 1, cex = 0.8), 
                              panel.scales)

    theta.par <- getPlotArgs("axis.line", local.resets = panel.scales,
                             user.resets = grid, elements = "theta",
                             defaults.only = FALSE)
   
##################
#these should be rationalised
##################

    if(!is.null(theta.par$theta))
        grid.theta <- theta.par$theta
    if(!is.null(theta.par$grid.theta))
        grid.theta <- theta.par$grid.theta
    if(is.null(grid.theta)){
        grid.theta <- if(!is.null(theta.par$at))
                           theta.par$at else
                               if(!is.null(thetalim))
                                   pretty(thetalim, theta.par$tick.number) else
                                      (0:12) * 30    
    }

    r.par <- getPlotArgs("axis.line", local.resets = panel.scales,
                         user.resets = grid, elements = "r", 
                         defaults.only=FALSE)

    
##need to look at this 
##maybe also above regarding the interval assignment for theta.
##little weird



    if(!is.null(r.par$r))
        grid.r <- r.par$r
    if(!is.null(r.par$grid.r))
        grid.r <- r.par$grid.r
    if(is.null(grid.r)){
        grid.r <- if(!is.null(r.par$at))
                        r.par$at else
                            if(!is.null(rlim))
                                pretty(rlim, r.par$tick.number) else
                                      NULL    
    }


    #make the arc grid.r 


    for(i in grid.r){
        theta <- seq(min(grid.theta), max(grid.theta), 
                     length.out = 10 * (max(grid.theta)-min(grid.theta)))
        r <- rep(i, length(theta))
        x <- r * sin(pi * theta/180)
        y <- r * cos(pi * theta/180)
        temp <- listUpdate(list(x=x, y=y, col="green"), theta.par)
        temp$col.line <- temp$col
        do.call(panel.lines, temp)
    }

##minor origin fix
##same in panel.polarAxes

    for(i in grid.theta){
        theta <- c(0, i)
#        r <- c(min(grid.r), max(grid.r))
        r <- c(0, max(grid.r)) 
        x <- r * sin(pi * theta/180)
        y <- r * cos(pi * theta/180)
        temp <- listUpdate(list(x=x, y=y, col="grey"), r.par)
        temp$col.line <- temp$col
        do.call(panel.lines, temp)
    }

}


#####################################################
#####################################################


panel.polarLabels <- function(labels.theta = NULL, labels.r = NULL,
         thetalim = NULL, rlim = NULL, ..., 
         labels = NULL, panel.scales = NULL){

    extra.args <- list(...)

    if (!is.list(panel.scales)) panel.scales <- list()
    if (!is.list(labels)) labels <- list()

    panel.scales <- listUpdate(list(draw = TRUE, arrows = FALSE, tick.number = 5, 
                                    abbreviate = FALSE, minlength = 4, tck = 1, 
                                    col.line = 1, cex = 0.8), 
                              panel.scales)


    theta.par <- getPlotArgs("axis.text", local.resets = panel.scales,
                             user.resets = labels, elements = "theta",
                             defaults.only = FALSE)

    
    if(!is.null(theta.par$theta))
        labels.theta <- theta.par$theta
    if(!is.null(theta.par$grid.theta))
        labels.theta <- theta.par$grid.theta
    if(is.null(labels.theta)){
        labels.theta <- if(!is.null(theta.par$at))
                           theta.par$at else
                               if(!is.null(thetalim))
                                   pretty(thetalim, theta.par$tick.number) else
                                      (0:12) * 30    
    }

    r.par <- getPlotArgs("axis.text", local.resets = panel.scales,
                         user.resets = labels, elements = "r", 
                         defaults.only=FALSE)



    if(!is.null(r.par$r))
        labels.r <- r.par$r
    if(!is.null(r.par$grid.r))
        labels.r <- r.par$grid.r
    if(is.null(labels.r)){
        labels.r <- if(!is.null(r.par$at))
                        r.par$at else
                            if(!is.null(rlim))
                                pretty(rlim, r.par$tick.number) else
                                      NULL    
    }

    #make the arc grid.r 

    for(i in labels.theta){

        theta <- rep(i, length(labels.r))
        x <- labels.r * sin(pi * theta/180)
        y <- labels.r * cos(pi * theta/180)
        temp <- listUpdate(list(x=x, y=y, labels=labels.r), theta.par)
        temp$col.line <- temp$col
        do.call(panel.text, temp)
    }

}
 


