#' Coefplot plotting
#'
#' Build ggplot object for coefplot
#'
#' This function builds up the ggplot layer by layer for \code{\link{coefplot.lm}}
#'
#' @author Jared P. Lander www.jaredlander.com
#' @seealso \code{\link{coefplot.default}} \code{\link{coefplot}} \code{\link{multiplot}}
#' @aliases buildPlotting.default
#' @param modelCI An object created by \code{\link{buildModelCI}}
#' @param title The name of the plot, if NULL then no name is given
#' @param xlab The x label
#' @param ylab The y label
#' @param innerCI How wide the inner confidence interval should be, normally 1 standard deviation.  If 0, then there will be no inner confidence interval.
#' @param outerCI How wide the outer confidence interval should be, normally 2 standard deviations.  If 0, then there will be no outer confidence interval.
#' @param multi logical; If this is for \code{\link{multiplot}} then leave the colors as determined by the legend, if FALSE then make all colors the same
#' @param lwdInner The thickness of the inner confidence interval
#' @param lwdOuter The thickness of the outer confidence interval
#' @param pointSize Size of coefficient point
#' @param color The color of the points and lines
#' @param shape The shape of the points
#' @param linetype The linetype of the error bars
#' @param cex The text size multiplier, currently not used
#' @param textAngle The angle for the coefficient labels, 0 is horizontal
#' @param numberAngle The angle for the value labels, 0 is horizontal
#' @param zeroColor The color of the line indicating 0
#' @param zeroLWD The thickness of the 0 line
#' @param zeroType The type of 0 line, 0 will mean no line
#' @param facet logical; If the coefficients should be faceted by the variables, numeric coefficients (including the intercept) will be one facet
#' @param scales The way the axes should be treated in a faceted plot.  Can be c("fixed", "free", "free_x", "free_y")
#' @param numeric logical; If true and factors has exactly one value, then it is displayed in a horizontal graph with constinuous confidence bounds.
#' @param fillColor The color of the confidence bounds for a numeric factor
#' @param alpha The transparency level of the numeric factor's confidence bound
#' @param horizontal logical; If the plot should be displayed horizontally
#' @param value Name of variable for value metric
#' @param coefficient Name of variable for coefficient names
#' @param errorHeight Height of error bars
#' @param dodgeHeight Amount of vertical dodging
#' @return a ggplot graph object
#' @examples
#'
#' data(diamonds)
#' model1 <- lm(price ~ carat + cut, data=diamonds)
#' theCI <- coefplot:::buildModelCI(model1)
#' coefplot:::buildPlotting.default(theCI)
#' coefplot(model1)
#'
buildPlotting.default <- function(modelCI, 
                                  title="Coefficient Plot", 
                                  xlab="Value", ylab="Coefficient", lwdInner=1, lwdOuter=0, pointSize=3,
                                  color="blue", cex=.8, textAngle=0, numberAngle=0, 
                                  shape=16, linetype=1,
                                  outerCI=2, innerCI=1, multi=FALSE, 
                                  zeroColor="grey", zeroLWD=1, zeroType=2, 
                                  numeric=FALSE, fillColor="grey", alpha=1/2,
                                  horizontal=FALSE, facet=FALSE, scales="free",
                                  value="Value", coefficient="Coefficient", 
                                  errorHeight=0, dodgeHeight=1)
{
    ## build the layer infos
    outerCIGeom <- geom_errorbarh(aes_string(xmin="LowOuter", xmax="HighOuter", color="Model", linetype="Model"), lwd=lwdOuter, height=errorHeight, position=position_dodgev(height=dodgeHeight))
    
    innerCIGeom <- geom_errorbarh(aes_string(xmin="LowInner", xmax="HighInner", color="Model", linetype="Model"),lwd=lwdInner, height=errorHeight, position=position_dodgev(height=dodgeHeight))
    
    # ribbon layer
    #ribbonGeom <- list(None=NULL, geom_ribbon(aes(ymin=LowOuter, ymax=HighOuter, group=Checkers), data=modelCI, fill=fillColor, alpha=alpha, lwd=lwdOuter))
    
    # point layer
    pointGeom <- geom_point(aes_string(xmin=value, xmax=value, color="Model", shape="Model"), size=pointSize, position=position_dodgev(height=dodgeHeight))

    #colorAes <- list(None=NULL, Single=aes(color=as.factor(Model)))
    colorScaleSingle <- scale_color_manual(values=rep(color, length(unique(modelCI$Model))), guide=FALSE)
    shapeScaleSingle <- scale_shape_manual(values=rep(shape, length(unique(modelCI$Model))), guide=FALSE)
    linetypeScaleSingle <- scale_linetype_manual(values=rep(linetype, length(unique(modelCI$Model))), guide=FALSE)
    
    xScale <- list(None=NULL, Single=scale_x_discrete())
    
    # faceting info
    faceting <- list(None=NULL, Display=facet_wrap(~Checkers, scales=scales))
    
    # for a regular coefplot or a multiplot in seperate facets
    #p <- ggplot(data=modelCI, aes_string(x=value))
    p <- ggplot(data=modelCI, aes_string(x=value, y=coefficient))    		# the basics of the plot
    #p <- p + colorAes[[1 + multi]] #                                    # in case model needs to be factorized, do it here
    p <- p + geom_vline(xintercept=0, colour=zeroColor, linetype=zeroType, lwd=zeroLWD)		# the zero line
    p <- p + outerCIGeom +    				# the outer CI bars
        innerCIGeom						# the inner CI bars
    p <- p + pointGeom						# the points
    #p <- p + xScale[[1 + multi]]
    p <- p + theme(axis.text.y=element_text(angle=textAngle, hjust=1), axis.text.x=element_text(angle=numberAngle, vjust=.5)) + 
        labs(title=title, x=xlab, y=ylab)    # labeling and text info
    p <- p + if(!multi){ list(colorScaleSingle, shapeScaleSingle, linetypeScaleSingle) }
    p <- p + faceting[[facet + 1]]    	# faceting
    p <- p + if(horizontal) coord_flip()
    
    return(p)		# return the ggplot object
}
