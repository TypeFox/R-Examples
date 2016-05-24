
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  .emptyPlot                Creates an empty plot page
#  .subPlot                  Creates a sub plot
#  .cnvrt.coords             Converts coordinates
################################################################################


.emptyPlot <- 
function()
{
    # Description:
    #   Creates an empty plot page
    
    # FUNCTION:
    
    # Plot:
    plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", 
        ann = TRUE, axes = FALSE)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
 

# Package: TeachingDemos
# Title: Demonstrations for teaching and learning
# Version: 1.6
# Author: Greg Snow
# Description: This package is a set of demonstration functions that can 
#   be used in a classroom to demonstrate statistical concepts, or on your 
#   own to better understand the concepts or the programming.
# Maintainer: Greg Snow <greg.snow@intermountainmail.org>
# License: Artistic-2.0
# Date: 2007-09-07
# Suggests: tcltk, lattice, tkrplot, MASS, rgl, sgeostat, mapproj
# Enhances: maptools, qcc
# LazyData: true
# Packaged: Wed Nov 28 12:13:53 2007; hornik


# ------------------------------------------------------------------------------


# Example:
#   plot(c(0, 1), c(0, 1)); .subPlot(hist(rnorm(100)), x = 0.5, y = 0.5)


# ------------------------------------------------------------------------------


.subPlot <- 
function(fun, x, y = NULL, size = c(1, 1), vadj = 0.5, hadj = 0.5, 
    pars = NULL)
{
    # Description:
    #   Creates a sub plot
    
    # FUNCTION:
    
    # Par:
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    
    if(missing(x)) x <- locator(2)
    xy <- xy.coords(x,y)
    
    if(length(xy$x) != 2){
    pin <- par('pin')
    tmp <- .cnvrt.coords(xy$x[1],xy$y[1],'usr')$plt
    
    x <- c( tmp$x - hadj*size[1]/pin[1],
            tmp$x + (1-hadj)*size[1]/pin[1] )
    y <- c( tmp$y - vadj*size[2]/pin[2],
            tmp$y + (1-vadj)*size[2]/pin[2] )
    
    xy <- .cnvrt.coords(x,y,'plt')$fig
    } else {
    xy <- .cnvrt.coords(xy,,'usr')$fig
    }
    
    par(pars)
    par(plt=c(xy$x,xy$y), new = TRUE)
    fun
    tmp.par <- par(no.readonly = TRUE)
    
    return(invisible(tmp.par))
}


# ------------------------------------------------------------------------------


.cnvrt.coords <- 
function(x,y=NULL,input=c('usr','plt','fig','dev','tdev')) 
{
    # Description:
    #   Converts Coordinates
    
    # FUNCTION:
    
    input <- match.arg(input)
    xy <- xy.coords(x,y)
  
    cusr <- par('usr')
    cplt <- par('plt')
    cfig <- par('fig')
    cdin <- par('din')
    comi <- par('omi')
    cdev <- c(comi[2]/cdin[1],(cdin[1]-comi[4])/cdin[1],
            comi[1]/cdin[2],(cdin[2]-comi[3])/cdin[2])
  
    if(input=='usr'){
        usr <- xy
        
        plt <- list()
        plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
        plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
        
        fig <- list()
        fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
        fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
        
        dev <- list()
        dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
        dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]
        
        tdev <- list()
        tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
        tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
        
        return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
    }

    if(input=='plt') {

        plt <- xy
        
        usr <- list()
        usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
        usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
        
        fig <- list()
        fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
        fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
        
        dev <- list()
        dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
        dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]
        
        tdev <- list()
        tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
        tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
        
        return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
    }

    if(input=='fig') {

        fig <- xy
        
        plt <- list()
        plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
        plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])
        
        usr <- list()
        usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
        usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
        
        dev <- list()
        dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
        dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]
        
        tdev <- list()
        tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
        tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
        
        return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
    }

    if(input=='dev'){
        dev <- xy
        
        fig <- list()
        fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
        fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])
        
        plt <- list()
        plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
        plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])
        
        usr <- list()
        usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
        usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
        
        tdev <- list()
        tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
        tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]  
        
        return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
    }

    if(input=='tdev'){
        tdev <- xy
    
        dev <- list()
        dev$x <- (tdev$x-cdev[1])/(cdev[2]-cdev[1])
        dev$y <- (tdev$y-cdev[3])/(cdev[4]-cdev[3])
        
        fig <- list()
        fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
        fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])
    
        plt <- list()
        plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
        plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])
    
        usr <- list()
        usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
        usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
    
        tdev <- list()
        tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
        tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]  
    
        return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
    }
  
}


################################################################################

