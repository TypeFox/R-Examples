
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:               PARAMETER ESTIMATION:
#  .gogarchFit             Fits the parameters of a GO-GARCH process
################################################################################


.gogarchFit <- 
function(formula = ~ garch(1, 1), data, 
    init.rec = c("mci", "uev"), 
    delta = 2, 
    skew = 1, 
    shape = 4, 
    cond.dist = c("norm", "snorm", "ged", "sged", "std", "sstd", "snig", "QMLE"), 
    include.mean = TRUE, 
    include.delta = NULL, 
    include.skew = NULL, 
    include.shape = NULL, 
    leverage = NULL, 
    trace = TRUE, 
    algorithm = c("nlminb", "lbfgsb", "nlminb+nm", "lbfgsb+nm"), 
    hessian = c("ropt", "rcd"), 
    control = list(), 
    title = NULL, 
    description = NULL, 
    ...) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Fits a GO-Garch Model using Independent Component Analysis
    
    # Arguments:
    #   The arguments are the same as for the univariate case.
    #   formula - formula for all marginal models
    #   data - multivariate timeSeries object
    #   ...
    
    # Value:
    #   S4 Object of class (univariate) fGARCH ...
    
    # Notes:
    #   This function has still a preliminary status ...
    #   This function was inspired from the contributed gogarch 
    #   package of Bernhard Pfaff.
    
    # Example:
    #   require(fEcofin); data(DowJones30)
    #   X = returns(as.timeSeries(DowJones30)); head(X)
    #   N = 5; ans = .gogarchFit(data = X[, 1:N], trace = FALSE); ans
    #   ans@h.t

    # FUNCTION:
      
    # Multivariate ?
    stopifnot(isMultivariate(data))
    
    
    # Data:
    X = data
    
    # Marginal Garch Models:
    garchControl = list(
        init.rec = init.rec, delta = delta, skew = skew, shape = shape, 
        cond.dist = cond.dist, include.mean = include.mean, 
        include.delta = include.delta, include.skew = include.skew, 
        include.shape = include.shape, leverage = leverage, 
        trace = trace, algorithm = algorithm, hessian = hessian, 
        control = control, title = title, description = description)
        
    # Compute fastICA:
    #   ... the following lines of code were borrowed from 
    #   Bernhard Pfaff's contributed package gogarch
    V <- t(X) %*% X / nrow(X)
    svd <- svd(V)
    P <- svd$u
    Dsqr <- diag(sqrt(svd$d))
    # set.seed(4711) 
    ica <- fastICA(X, n.comp = ncol(X))
    Z <- P %*% Dsqr %*% t(P) %*% ica$W
    colnames(Z) = rownames(Z) = colnames(data)
    Y <- X %*% solve(Z)

    # Fit Marginal Garch Models:
    fit <- apply(Y, 2, function(x) do.call("garchFit", 
        c(list(formula = formula, data = x), garchControl)))
    
    # Compute Conditional Variances:
    #   ... the following lines of code were borrowed from 
    #   Bernhard Pfaff's contributed package gogarch
    H <- matrix(unlist(lapply(fit, function(x) x@h.t)), 
        ncol = ncol(X), nrow = nrow(X))
    Hdf = data.frame(t(H))
    rownames(Hdf) = colnames(data)
    colnames(Hdf) = rownames(data)
    H.t <- lapply(Hdf, function(x) Z %*% diag(x) %*% t(Z))
    
    # Add Title and Description:
    if(is.null(title)) title = "ICA GO-GARCH Modelling"
    if(is.null(description)) description = description()
    
    # Result:
    ans = new("fGARCH",
        call = as.call(match.call()),
        formula = formula,
        method = "ICA go-Garch Parmeter Estimation",
        data = c(Records = nrow(data), Instruments = ncol(data)),               
        fit = fit,
        residuals = numeric(),      
        fitted = numeric(),         
        h.t = c(Records = length(H.t), Dimension =dim(H.t[[1]])),             
        sigma.t = numeric(),         
        title = title,               
        description = description
    )
    
    # Multivariate Series:
    attr(ans@data, "data") <- data
    attr(ans@h.t, "H.t") <- H.t
    
    # Return Value:
    ans   
}


################################################################################

