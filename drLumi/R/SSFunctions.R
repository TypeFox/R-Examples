#' Self-Starting Nls exponential regression model
#' 
#' This selfStart model evaluates the exponential growth 
#' regression model and its gradient. It has an \code{initial}
#' attribute that will evaluate initial estimates of the parameters 
#' \code{y0}, and \code{b} for a given set of data. Instead of the standard 
#' \code{exp} function this implementation use the \code{10^} function.
#' \deqn{f(x)=y_0 \times 10^b}
#' @param x a numeric vector of values at which to evaluate the model
#' @param y0 a numeric parameter representing the value of the response when 
#' x is 0
#' @param b a numeric parameter representing the growth rate
#' @return The value returned is a list containing the nonlinear function, 
#' the self starter function and the parameter names.
#' @usage SSexp(x, b, y0)
#' 
#' @format A selfStart model
#' 
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' # Select analyte FGF for plate 1
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)$plate_1
#'
#' ig <- scluminex("plate_1",sdf$standard, sdf$background, 
#'             lfct="SSexp", 
#'             bkg="ignore", 
#'             fmfi="mfi", 
#'             verbose=FALSE)
#'
#' summary(ig)
#'
#' 
#' @author John Aponte <john.aponte@@cresib.cat>
#' @export
SSexp <- selfStart( ~ y0 * 10^(x/b),
            function(mCall, LHS, data)
            {
            xy <- sortedXyData(mCall[["x"]], LHS, data)
            xypos <- xy[which(xy[,"y"]>0),]
            lmFit <- lm(log10(xypos[,"y"])~xypos[,"x"])
            coefs <- coef(lmFit)
            y0 <- 10^(coefs[1])
            b <- 1/coefs[2]
            value <- c(b,y0)
            names(value) <- mCall[c("b","y0")]
            value
            }, c("b", "y0"))


# Inverse function for the SSexp  (x ~ fx)
# 
# @param fx MFI value
# @param b coefficient of b
# @param c coefficient of y0
inSSexp <- function(fx, b, c){
    ans <- list()  
    ans$est <- b * log(fx/c)/log(10)
    ans$formtext <- "~  x1 * log(fx/x2)/log(10)"
    ans$formtext.cons <-   "~  x1 * log(fx/cons)/log(10)"
    ans
} 

# First order derivative of the SSexp function
# 
# @param x concentration value
# @param b coefficient of b
# @param y0 coefficient of y0
dSSexp <- function(x, b, y0){
    ans <- y0 * (10^(x/b) * (log(10) * (1/b)))
    ans  
}

# Second order derivative of the SSexp function
#  
# @param x concentration value
# @param b coefficient of b
# @param y0 coefficient of y0
d2SSexp <- function(x, b, y0){
    ans <- y0 * (10^(x/b) * (log(10) * (1/b)) * (log(10) * (1/b)))
    ans  
}

# Third order derivative of the SSexp function
# 
# @param x concentration value
# @param b coefficient of b
# @param y0 coefficient of y0
d3SSexp <- function(x, b, y0){
    ans <- y0 * (10^(x/b) * (log(10) * (1/b)) * (log(10) * (1/b)) 
            * (log(10) * (1/b)))
    ans  
}


#' Self-Starting Nls 4 parameters logistic regression model
#' 
#' This selfStart model evaluates the 4 parameters logistic  
#' regression model and its gradient.  It has an \code{initial}
#' attribute that will evaluate initial estimates of the parameters 
#' \code{hAsym}, \code{lAsym}, \code{Slope} and \code{xMid} 
#' for a given set of data
#' Instead of the standard \code{exp} function this implementation use 
#' the \code{10^} function.
#' \deqn{f(x)=lAsym +\frac{hAsym-lAsym}{1+10^{Slope(x-xMid)}}}
#' @param x a numeric vector of values at which to evaluate the model
#' @param hAsym a numeric parameter representing the higher asymptote when 
#' \code{x} trend to \code{Inf}
#' @param lAsym a numeric parameter representing the lower asymptote when 
#' \code{x} trend to \code{-Inf}
#' @param xMid is the x value corresponding to the inflection point
#' @param Slope a numeric parameter representing the \code{-slope} of the 
#' function at the inflection point
#' @return The value returned is a list containing the nonlinear function, 
#' the self starter function and the parameter names.
#' @usage SSl4(x, Slope, lAsym, hAsym, xMid)
#' @format A selfStart model
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' # Select analyte FGF for plate 1
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)[[1]]
#'
#' ig <- scluminex("plate_1",sdf$standard, sdf$background, 
#'            lfct="SSl4", 
#'            bkg="ignore", 
#'            fmfi="mfi", 
#'            verbose=FALSE)
#'  
#' summary(ig)
#' 
#' 
#' @author John Aponte <john.aponte@@cresib.cat>
#' @export
SSl4 <- selfStart( ~ lAsym + (hAsym-lAsym)/(1+10^(Slope*(x-xMid))),
            function(mCall, LHS, data)
            {
            xy <- sortedXyData(mCall[["x"]], LHS, data)
            plowl <- NLSstLfAsymptote(xy)
            phighl <- NLSstRtAsymptote(xy)
            pslope <- -1
            ped50 <- NLSstClosestX(xy, median(xy[,"y"]))
            value <- c(pslope,plowl,phighl,ped50)
            names(value)<-mCall[c("Slope","lAsym","hAsym","xMid")]
            value
            },c("Slope","lAsym","hAsym","xMid"))




# Inverse function for the SSl4  (x ~ fx)
# 
# @param fx MFI value
# @param b slope parameter
# @param c low asymptote parameter
# @param d high asymptote parameter
# @param e mid concentration parameter
inSSl4 <- function(fx, b, c, d, e){
    ans <- list()
    ans$est <- (log((d - fx )/(fx - c)) * 1/(log(10)*b)) + e
    ans$formtext <- "~ (log((x3 - fx )/(fx - x2)) * 1/(log(10)*x1)) + x4"
    ans$formtext.cons <- "~ (log((x2 - fx )/(fx - cons)) 
                            * 1/(log(10)*x1)) + x3"
    ans
}

# First order derivative of the SSl4 function
# 
# @param x a numeric vector of values at which to evaluate the model
# @param Slope a numeric parameter representing the \code{-slope} of the 
# function at the inflection point
# @param lAsym a numeric parameter representing the lower asymptote when 
# \code{x} trend to \code{-Inf}
# @param hAsym a numeric parameter representing the higher asymptote when 
# \code{x} trend to \code{Inf}
# @param xMid is the x value corresponding to the inflection point
dSSl4 <- function(x, Slope, lAsym, hAsym, xMid){
    ans <- -((hAsym - lAsym) * (10^(Slope * (x - xMid)) * 
            (log(10) * Slope))/
            (1 + 10^(Slope * (x - xMid)))^2)
    ans
}

# Second order derivative of the SSl4 function
# 
# @param x a numeric vector of values at which to evaluate the model
# @param Slope a numeric parameter representing the \code{-slope} of the 
# function at the inflection point
# @param lAsym a numeric parameter representing the lower asymptote when 
# \code{x} trend to \code{-Inf}
# @param hAsym a numeric parameter representing the higher asymptote when 
# \code{x} trend to \code{Inf}
# @param xMid is the x value corresponding to the inflection point
d2SSl4 <- function(x, Slope, lAsym, hAsym, xMid){ 
    ans <- -((hAsym - lAsym) * (10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (log(10) * Slope))/(1 + 10^(Slope * (x - xMid)))^2 - (hAsym - 
            lAsym) * (10^(Slope * (x - xMid)) * (log(10) * Slope)) * (2 * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (1 + 10^(Slope * (x - xMid)))))/((1 + 10^(Slope * 
            (x - xMid)))^2)^2)
    ans
}

# Third order derivative of the SSl4 function
# 
# @param x a numeric vector of values at which to evaluate the model
# @param Slope a numeric parameter representing the \code{-slope} of the 
# function at the inflection point
# @param lAsym a numeric parameter representing the lower asymptote when 
# \code{x} trend to \code{-Inf}
# @param hAsym a numeric parameter representing the higher asymptote when 
# \code{x} trend to \code{Inf}
# @param xMid is the x value corresponding to the inflection point
d3SSl4 <- function(x, Slope, lAsym, hAsym, xMid){ 
    ans <- -((hAsym - lAsym) * (10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (log(10) * Slope) * (log(10) * Slope))/(1 + 10^(Slope * 
            (x - xMid)))^2 - 
            (hAsym - lAsym) * (10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (log(10) * Slope)) * (2 * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope) * (1 + 10^(Slope * 
            (x - xMid)))))/((1 + 
            10^(Slope * (x - xMid)))^2)^2 - (((hAsym - lAsym) * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope) * (log(10) * Slope)) * (2 * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope) * (1 + 10^(Slope * 
            (x - xMid))))) + (hAsym - lAsym) * (10^(Slope * (x - 
            xMid)) * (log(10) * Slope)) * (2 * (10^(Slope * (x - xMid)) * 
            (log(10) * Slope) * (log(10) * Slope) * (1 + 10^(Slope * 
            (x - xMid))) + 10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope)))))/((1 + 10^(Slope * 
            (x - xMid)))^2)^2 - (hAsym - lAsym) * (10^(Slope * (x - xMid)) * 
            (log(10) * Slope)) * (2 * (10^(Slope * (x - xMid)) * (log(10) * 
            Slope) * (1 + 10^(Slope * (x - xMid))))) * (2 * (2 * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope) * (1 + 10^(Slope * (x - xMid)))) * 
            ((1 + 10^(Slope * (x - xMid)))^2)))/(((1 + 10^(Slope * (x - 
            xMid)))^2)^2)^2))
    ans
}


#' Self-Starting Nls 5 parameters logistic regression model
#' 
#' This selfStart model evaluates the 5 parameters logistic  
#' regression model and its gradient. 
#' It has an \code{initial}
#' attribute that will evaluate initial estimates of the parameters 
#' \code{hAsym}, \code{lAsym}, \code{Slope}, \code{xMid} and \code{Asymetry} 
#' for a given set of data
#' Instead of the standard \code{exp} function this implementation use 
#' the \code{10^} function.
#' \deqn{f(x)=lAsym +\frac{hAsym-lAsym}{(1+10^{Slope(x-xMid)})^{Asymetry}}}
#' @param x a numeric vector of values at which to evaluate the model
#' @param hAsym a numeric parameter representing the higher asymptote 
#' when \code{x} trend to \code{Inf}
#' @param lAsym a numeric parameter representing the lower asymptote when 
#' \code{x} trend to \code{-Inf}
#' @param xMid is the x value corresponding to the inflection point
#' @param Slope is a numeric parameter representing the \code{-slope} of the 
#' function at the inflection point
#' @param Asymetry is a numeric parameter representing the asymetry around 
#' the inflection point
#' @return The value returned is a list containing the nonlinear function, 
#' the self starter function and the parameter names.
#' @usage SSl5(x, Slope, lAsym, hAsym, xMid, Asymetry)
#' @format A selfStart model
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' # Select analyte FGF for plate 1
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)[[1]]
#'
#' # SSl5 
#' ig <- scluminex("plate_1",sdf$standard, sdf$background, 
#'             lfct="SSl5", 
#'             bkg="ignore", 
#'             fmfi="mfi", 
#'             verbose=FALSE)
#' 
#' summary(ig)
#' 
#' @author John Aponte <john.aponte@@cresib.cat>
#' @export
SSl5 <- selfStart( ~ lAsym + (hAsym-lAsym)/(1+10^(Slope*(x-xMid)))^Asymetry,
        function(mCall, LHS, data)
        {
        xy <- sortedXyData(mCall[["x"]], LHS, data)
        plowl<- NLSstLfAsymptote(xy)
        phighl<-NLSstRtAsymptote(xy)
        pslope<- -1
        pf<-1
        ped50<-NLSstClosestX(xy, mean(xy[,"y"]))
        value<-c(pslope,plowl,phighl,ped50,pf)
        names(value)<-mCall[c("Slope","lAsym","hAsym","xMid","Asymetry")]
        value
        },c("Slope","lAsym","hAsym","xMid","Asymetry"))


# Inverse function for the SSl5  (x ~ fx)
# 
# @param fx MFI value
# @param b slope parameter
# @param c low asymptote parameter
# @param d high asymptote parameter
# @param e mid concentration parameter
# @param f asymetry parameter
inSSl5 <- function(fx, b, c, d, e, f){
    ans <- list()  
    ans$est <- (log( ((d-c)/(fx-c))^(1/f) - 1  ) * 1/(log(10)*b)) + e
    t1 <- "~ (log( ((x3-x2)/(fx-x2))^(1/x5) - 1  ) * 1/(log(10)*x1)) + x4"
    ans$formtext <- t1
    t2 <- "~ (log( ((x2-cons)/(fx-cons))^(1/x4) - 1  ) * 1/(log(10)*x1)) + x3"
    ans$formtext.cons <- t2
    ans
}

# First order derivative of the SSl5 function
# 
# @param x a numeric vector of values at which to evaluate the model
# @param Slope is a numeric parameter representing the \code{-slope} of the 
# function at the inflection point
# @param lAsym a numeric parameter representing the lower asymptote when 
# \code{x} trend to \code{-Inf}
# @param hAsym a numeric parameter representing the higher asymptote when 
# \code{x} trend to \code{Inf}
# @param xMid is the x value corresponding to the inflection point
# @param Asymetry is a numeric parameter representing the asymetry around 
# the inflection point
dSSl5 <- function(x, Slope, lAsym, hAsym, xMid, Asymetry){    
    ans <- -((hAsym - lAsym) * ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 
            1) * (Asymetry * (10^(Slope * (x - xMid)) * 
            (log(10) * Slope))))/((1 + 10^(Slope * (x - xMid)))^Asymetry)^2)
    ans
}

# Second order derivative of the SSl5 function
# 
# @param x a numeric vector of values at which to evaluate the model
# @param Slope is a numeric parameter representing the \code{-slope} 
# of the function at the inflection point
# @param lAsym a numeric parameter representing the lower asymptote 
# when \code{x} trend to \code{-Inf}
# @param hAsym a numeric parameter representing the higher asymptote 
# when \code{x} trend to \code{Inf}
# @param xMid is the x value corresponding to the inflection point
# @param Asymetry is a numeric parameter representing the asymetry 
# around the inflection point
d2SSl5 <- function(x, Slope, lAsym, hAsym, xMid, Asymetry){ 
    ans <-  -((hAsym - lAsym) * ((1 + 10^(Slope * (x - xMid)))^((Asymetry - 
            1) - 1) * ((Asymetry - 1) * (10^(Slope * (x - xMid)) * (log(10) * 
            Slope))) * (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * 
            Slope))) + (1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * 
            (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (log(10) * Slope))))/((1 + 10^(Slope * (x - xMid)))^Asymetry)^2 - 
            (hAsym - lAsym) * ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 
            1) * (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * 
            Slope)))) * (2 * ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 
            1) * (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * 
            Slope))) * ((1 + 10^(Slope * (x - xMid)))^Asymetry)))/(((1 + 
            10^(Slope * (x - xMid)))^Asymetry)^2)^2)
    ans
}

# Third order derivative of the SSl5 function
# 
# @param x a numeric vector of values at which to evaluate the model
# @param Slope is a numeric parameter representing the \code{-slope} of the 
# function at the inflection point
# @param lAsym a numeric parameter representing the lower asymptote when 
# \code{x} trend to \code{-Inf}
# @param hAsym a numeric parameter representing the higher asymptote when 
# \code{x} trend to \code{Inf}
# @param xMid is the x value corresponding to the inflection point
# @param Asymetry is a numeric parameter representing the asymetry around the 
# inflection point
d3SSl5 <- function(x, Slope, lAsym, hAsym, xMid, Asymetry){   
    ans <-  -((hAsym - lAsym) * (((1 + 10^(Slope * (x - xMid)))^(((Asymetry -
            1) - 1) - 1) * (((Asymetry - 1) - 1) * (10^(Slope * (x - 
            xMid)) * (log(10) * Slope))) * ((Asymetry - 1) * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope))) + (1 + 10^(Slope * (x - 
            xMid)))^((Asymetry - 1) - 1) * ((Asymetry - 1) * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope) * (log(10) * Slope)))) * 
            (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * Slope))) + 
            (1 + 10^(Slope * (x - xMid)))^((Asymetry - 1) - 1) * ((Asymetry - 
            1) * (10^(Slope * (x - xMid)) * (log(10) * Slope))) * 
            (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (log(10) * Slope))) + ((1 + 10^(Slope * (x - xMid)))^((Asymetry - 
            1) - 1) * ((Asymetry - 1) * (10^(Slope * (x - xMid)) * (log(10) * 
            Slope))) * (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * 
            Slope) * (log(10) * Slope))) + (1 + 10^(Slope * 
            (x - xMid)))^(Asymetry - 
            1) * (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * Slope) * 
            (log(10) * Slope) * (log(10) * Slope)))))/((1 + 10^(Slope * 
            (x - xMid)))^Asymetry)^2 - (hAsym - lAsym) * ((1 + 10^(Slope * 
            (x - xMid)))^((Asymetry - 1) - 1) * ((Asymetry - 1) * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope))) * (Asymetry * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope))) + (1 + 10^(Slope * (x - 
            xMid)))^(Asymetry - 1) * (Asymetry * (10^(Slope * (x - xMid)) * 
            (log(10) * Slope) * (log(10) * Slope)))) * (2 * ((1 + 10^(Slope * 
            (x - xMid)))^(Asymetry - 1) * (Asymetry * (10^(Slope * (x - 
            xMid)) * (log(10) * Slope))) * 
            ((1 + 10^(Slope * (x - xMid)))^Asymetry)))/(((1 + 
            10^(Slope * (x - xMid)))^Asymetry)^2)^2 - (((hAsym - lAsym) * 
            ((1 + 10^(Slope * (x - xMid)))^((Asymetry - 1) - 1) * ((Asymetry - 
            1) * (10^(Slope * (x - xMid)) * (log(10) * Slope))) * 
            (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * Slope))) + 
            (1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope) * (log(10) * 
            Slope)))) * (2 * ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 
            1) * (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * Slope))) * 
            ((1 + 10^(Slope * (x - xMid)))^Asymetry))) + (hAsym - lAsym) * 
            ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope)))) * (2 * 
            (((1 + 10^(Slope * (x - xMid)))^((Asymetry - 1) - 1) * ((Asymetry - 
            1) * (10^(Slope * (x - xMid)) * (log(10) * Slope))) * 
            (Asymetry * (10^(Slope * (x - xMid)) * (log(10) * Slope))) + 
            (1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope) * (log(10) * 
            Slope)))) * ((1 + 10^(Slope * (x - xMid)))^Asymetry) + 
            (1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope))) * 
            ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope)))))))/(((1 + 
            10^(Slope * (x - xMid)))^Asymetry)^2)^2 - (hAsym - lAsym) * 
            ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope)))) * (2 * 
            ((1 + 10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * 
            (10^(Slope * (x - xMid)) * (log(10) * Slope))) * ((1 + 
            10^(Slope * (x - xMid)))^Asymetry))) * (2 * (2 * ((1 + 
            10^(Slope * (x - xMid)))^(Asymetry - 1) * (Asymetry * (10^(Slope * 
            (x - xMid)) * (log(10) * Slope))) * ((1 + 10^(Slope * (x - 
            xMid)))^Asymetry)) * (((1 + 10^(Slope * 
            (x - xMid)))^Asymetry)^2)))/((((1 + 
            10^(Slope * (x - xMid)))^Asymetry)^2)^2)^2))
    ans
}





#' Self-Starting Nls exponential constraint regression model
#' 
#' This selfStart model evaluates the exponential growth 
#' regression model and its gradient. It has an \code{initial}
#' attribute that will evaluate initial estimates of the parameters 
#' \code{y0}, and \code{b} for a given set of data. Instead of the standard 
#' \code{exp} function this implementation use the \code{10^} function.
#' \deqn{f(x)=y_0 \times 10^b}
#' 
#'  
#' @param ..constraint.value a numeric value representing the value of the 
#' response when \code{x} is 0.  In this 
#' function this value is not a parameter is just a numeric value to 
#' constraint \code{y_0} parameter.  
#' @param x a numeric vector of values at which to evaluate the model
#' @param b a numeric parameter representing the growth rate
#' @return The value returned is a list containing the nonlinear function, 
#' the self starter function and the parameter names.
#' 
#' @usage SSexpcons(..constraint.value, x, b)
#' 
#' @format A selfStart model
#' 
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' # Select analyte FGF for plate 1
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)$plate_1
#'
#' cons <- scluminex("plate_1",sdf$standard, sdf$background, 
#'             lfct="SSexp", 
#'             bkg="constraint", 
#'             fmfi="mfi", 
#'             verbose=FALSE)
#'
#' summary(cons)
#' 
#' # Comparison constraint vs no constraint (same returning value but estimate
#' # one parameter).
#' b <- 3
#' y0 <- 1
#' concentration <- 2
#' SSexp(concentration, b, y0)
#' SSexpcons(y0, concentration, b)
#'
#' @export
SSexpcons <- selfStart( ~ ..constraint.value * 10^(x/b),
                    function(mCall, LHS, data)
                    {
                      xy <- sortedXyData(mCall[["x"]], LHS, data)
                      xypos <- xy[which(xy[,"y"]>0),]
                      lmFit <- lm(log10(xypos[,"y"])~xypos[,"x"])
                      coefs <- coef(lmFit)
                       b <- 1/coefs[2]
                      value <- c(b)
                      names(value) <- mCall[c("b")]
                      value
                    }, c("b"))


#' Self-Starting Nls 4 parameters logistic constraint regression model
#' 
#' This selfStart model evaluates the 4 parameters logistic  
#' regression model and its gradient.  It has an \code{initial}
#' attribute that will evaluate initial estimates of the parameters 
#' \code{hAsym}, \code{Slope} and \code{xMid} 
#' for a given set of data.
#' Instead of the standard \code{exp} function this implementation use 
#' the \code{10^} function.
#' \deqn{f(x)=lAsym +\frac{hAsym-lAsym}{1+10^{Slope(x-xMid)}}}
#' 
#' 
#' 
#' @param ..constraint.value a numeric value representing 
#' the lower asymptote when  \code{x} trend to \code{-Inf}. In this 
#' function this value is not a parameter is just a numeric value to 
#' constraint \code{lAsym} parameter. 
#' @param x a numeric vector of values at which to evaluate the model
#' @param hAsym a numeric parameter representing the higher asymptote when 
#' \code{x} trend to \code{Inf}
#' @param xMid is the x value corresponding to the inflection point
#' @param Slope a numeric parameter representing the \code{-slope} of the 
#' function at the inflection point
#' 
#' 
#' @return The value returned is a list containing the nonlinear function, 
#' the self starter function and the parameter names.
#' 
#' @usage SSl4cons(..constraint.value, x, Slope, hAsym, xMid)
#' @format A selfStart model
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' # Select analyte FGF for plate 1
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)[[1]]
#'
#' cons <- scluminex("plate_1",sdf$standard, sdf$background, 
#'            lfct="SSl4", 
#'            bkg="constraint", 
#'            fmfi="mfi", 
#'            verbose=FALSE)
#'  
#' summary(cons)
#' 
#' # Comparison constraint vs no constraint (same returning value but 
#' # estimate 3 parameters).
#' lAsym <- 1
#' Slope <- 2
#' hAsym   <- 2
#' xMid <- 3
#' concentration <- 2
#' SSl4(concentration, Slope, lAsym, hAsym, xMid)
#' SSl4cons(lAsym, concentration, Slope, hAsym, xMid)
#' 
#' @export
SSl4cons <- selfStart( ~ ..constraint.value + (hAsym-..constraint.value)/(1+10^(Slope*(x-xMid))),
                       function(mCall, LHS, data)
                       {
                         xy <- sortedXyData(mCall[["x"]], LHS, data)
                         phighl <- NLSstRtAsymptote(xy)
                         pslope <- -1
                         ped50 <- NLSstClosestX(xy, median(xy[,"y"]))
                         value <- c(pslope,phighl,ped50)
                         names(value)<-mCall[c("Slope","hAsym","xMid")]
                         value
                       },c("Slope","hAsym","xMid"))



#' Self-Starting Nls 5 parameters logistic constraint regression model 
#' 
#' 
#' This selfStart model evaluates the 5 parameters logistic  
#' regression model and its gradient for the lower asymptote constraint method. 
#' It has an \code{initial}
#' attribute that will evaluate initial estimates of the parameters 
#' \code{hAsym}, \code{Slope}, \code{xMid} and \code{Asymetry} 
#' for a given set of data
#' Instead of the standard \code{exp} function this implementation use 
#' the \code{10^} function.
#' \deqn{f(x)=lAsym +\frac{hAsym-lAsym}{(1+10^{Slope(x-xMid)})^{Asymetry}}}
#' 
#'
#' @param ..constraint.value a numeric value representing 
#' the lower asymptote when  \code{x} trend to \code{-Inf}. In this 
#' function this value is not a parameter is just a numeric value to 
#' constraint \code{lAsym} parameter. 
#' @param x a numeric vector of values at which to evaluate the model
#' @param hAsym a numeric parameter representing the higher asymptote 
#' when \code{x} trend to \code{Inf}
#' @param xMid is the x value corresponding to the inflection point
#' @param Slope is a numeric parameter representing the \code{-slope} of the 
#' function at the inflection point
#' @param Asymetry is a numeric parameter representing the asymetry around 
#' the inflection point
#' 
#' 
#' @return The value returned is a list containing the nonlinear function, 
#' the self starter function and the parameter names.
#' 
#' 
#' @usage SSl5cons(..constraint.value,x, Slope, hAsym, xMid, Asymetry)
#' @format A selfStart model
#' 
#' 
#' @examples
#' # Load data
#' data(ecdata)
#' data(mfidata)
#' 
#' # Select analyte FGF for plate 1
#' dat <- mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",]
#' 
#' sdf <- data_selection(dat, ecdata)[[1]]
#'
#' # SSl5 
#' cons <- scluminex("plate_1",sdf$standard, sdf$background, 
#'             lfct="SSl5", 
#'             bkg="constraint", 
#'             fmfi="mfi", 
#'             verbose=FALSE)
#' 
#' summary(cons)
#' 
#' # Comparison constraint vs no constraint (same returning value but estimate
#' # 4 parameters).
#' lAsym <- 1
#' Slope <- 2
#' hAsym   <- 2
#' xMid <- 3
#' Asymetry <- 1.5
#' 
#' concentration <- 2
#' SSl5(concentration, Slope, lAsym, hAsym, xMid, Asymetry)
#' SSl5cons(lAsym, concentration, Slope, hAsym, xMid, Asymetry)
#' 
#' @export
SSl5cons <- selfStart( ~ ..constraint.value + (hAsym-..constraint.value)/(1+10^(Slope*(x-xMid)))^Asymetry,
                   function(mCall, LHS, data)
                   {
                     xy <- sortedXyData(mCall[["x"]], LHS, data)
                     phighl<-NLSstRtAsymptote(xy)
                     pslope<- -1
                     pf<-1
                     ped50<-NLSstClosestX(xy, mean(xy[,"y"]))
                     value<-c(pslope,phighl,ped50,pf)
                     names(value)<-mCall[c("Slope","hAsym","xMid","Asymetry")]
                     value
                   },c("Slope","hAsym","xMid","Asymetry"))



