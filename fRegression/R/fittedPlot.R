
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA. 


################################################################################
# FUNCTION:                   REGRESSION TERM PLOTS:
#  .fittedPlot                 Line Plot          
#  .fittedPersp                Perspective Plot         
#  .fittedContour              Contour Plot             
################################################################################


.fittedPlot <- 
    function(object, which = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    
    # Arguments:
    #   object - an object of class 'fREG' as returned by the function regFit
    
    # FUNCTION:
    
    model = object@fit$model
    responseName = colnames(model)[attr(terms(object), "response")]
    model.mat = as.matrix(object@fit$model)[,-attr(terms(object), "response")]
    N = NCOL(model.mat)
    zero = rep(0, times = N)
    
    if (is.null(which)) which = 1:N
    colNames = colnames(model.mat)[which]
    
    ans = NULL
    for (i in which) {
        one = zero
        one[i] = 1
        new.model.mat = model.mat
        new.model.mat = 0 * model.mat
        x = new.model.mat[, i] = model.mat %*% one
        y = predict(object, newdata = as.data.frame(new.model.mat))
        ans = cbind(ans, y)
        plot(x, y, xlab = colNames[i], ylab = paste("Fitted", colNames[i]))
    }
    
    colnames(ans) = paste(responseName, "(", colNames, ")", sep = "")
    as.data.frame(ans)
}


# ------------------------------------------------------------------------------   


.fittedPersp <- 
    function(object)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    
    # Arguments:
    #   object - an object of class 'fREG' as returned by the function regFit
    
    # FUNCTION:
    
    # Settings:
    model = object@fit$model
    responseName = colnames(model)[attr(terms(object), "response")]
    model.mat = as.matrix(object@fit$model)[,-attr(terms(object), "response")]
    N = NCOL(model.mat)
    colNames = colnames(model.mat)
    
    for (i in 1:(N-1)) {    
        rangeX = range(model.mat[, i])   
        X = seq(rangeX[1], rangeX[2], length = 10)
        newdata = matrix(rep(0, times = N*10*10), ncol = N)
        newdata[ ,i] = X
        for (j in (i+1):N) { 
            rangeY = range(model.mat[, j])   
            Y = seq(rangeY[1], rangeY[2], length = 10)
            XY = gridVector(X, Y)
            newdata[, j] = Y
            colnames(newdata) = colNames
            print(head(newdata))
            
            Z = predict(object, as.data.frame(newdata))
            Z = matrix(Z, ncol = 10)
            .perspPlot(X, Y, Z, xlab = colNames[i], ylab = colNames[j])
        }
    }

}


# ------------------------------------------------------------------------------


.fittedContour <- 
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    
    # Arguments:
    #   object - an object of class 'fREG' as returned by the function regFit
    
    # FUNCTION:
    
    # Settings:
    model <- object@fit$model
    responseName = colnames(model)[attr(terms(object), "response")]
    model.mat = as.matrix(object@fit$model)[,-attr(terms(object), "response")]
    N = NCOL(model.mat)
    colNames = colnames(model.mat)
    
    for (i in 1:(N-1)) {    
        rangeX = range(model.mat[, i])   
        X = seq(rangeX[1], rangeX[2], length = 10)
        newdata = matrix(rep(0, times = N*10*10), ncol = N)
        newdata[ ,i] = X
        for (j in (i+1):N) { 
            rangeY = range(model.mat[, j])   
            Y = seq(rangeY[1], rangeY[2], length = 10)
            XY = gridVector(X, Y)
            newdata[, j] = Y
            colnames(newdata) = colNames
            print(head(newdata))
            
            Z = predict(object, as.data.frame(newdata))
            Z = matrix(Z, ncol = 10)
            .contourPlot(X, Y, Z, xlab = colNames[i], ylab = colNames[j])
        }
    }

}


################################################################################    


