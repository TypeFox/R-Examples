
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

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
#  FUNCTION:               DESCRIPTION:
#
#  gsSelect                Find best fitted model according to AIC criterion
################################################################################

gsSelect <-
  function(
    data,
    order.max = c(1,1,1,1),
    selection.criteria = c("AIC", "AICc", "BIC"),
    is.aparch = FALSE,
    cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
    include.mean = TRUE, 
    algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"),
    ...)
{
      
    
    # TESTING
    #data = x; order.max = c(1,1,1,1); is.aparch = TRUE; cond.dist = "norm"; algorithm = "sqp"
    # Description:
    #   Estimate several models using GSgarch.Fit function and evaluates
    #   the AIC to decide which one is the best model
    
    # Arguments:
    #   data - vector of data
    #   intercept - a logical, should the mean value be estimated ? 
    #   nMAX, mMAX, pMAX, qMAX - maximum order to be estimated
    #   for the model ARMA(m,n)-GARCH/APARCH(p,q)
    #   cond.dist - name of the conditional distribution, one of
    #       gev, stable, norm, std, sstd
    #   algorithm - 
    #   APARCH - boolean, indicates whether the model is an APARCH model or not
    
    # Return:
    #   fit.min - The best model return from function GSgarch.Fit
    
    # FUNCTION:      
      
    # error treatment on input parameters
    cond.dist = match.arg(cond.dist)
    algorithm = match.arg(algorithm)
    selection.criteria = match.arg(selection.criteria)
    if( !is.numeric(data) || !is.vector(data))
        stop("data set must be a numerical one dimensional vector")
    
    # begin of function
    goodness.of.fit.min <- 1e99
    goodness.of.fit.current <- 1e99
    fit.min <- list()
    fit <- list()
    order.list <- .getOrder(order.max = order.max)
    order.list.size <- length(order.list[,1])
    for( i in 1:order.list.size)
    {
        m = order.list[i,1]; n = order.list[i,2]; p = order.list[i,3]; q = order.list[i,4]
        
        if(is.aparch == TRUE)
        {   
            formula = as.formula(paste ("~ arma(",m,", ",n,") + aparch(",p,", ",q,")", sep = "", 
                                      collapse = NULL))
        }
        else {   
          formula = as.formula(paste ("~ arma(",m,", ",n,") + garch(",p,", ",q,")", sep = "", 
                                      collapse = NULL))
        }
        cat("\n------------------------------------------------------------------------------------------\n")
        cat(paste("Model: ", paste(formula)[2], " with '", 
            cond.dist, "' conditional distribution",sep = ""))
        cat("\n------------------------------------------------------------------------------------------\n")
        fit = gsFit (data = data, formula = formula, cond.dist = cond.dist, 
                     algorithm = algorithm, include.mean = include.mean, ...)
        
        goodness.of.fit.current = fit@fit$ics[selection.criteria]
        if (goodness.of.fit.current < goodness.of.fit.min)
        {
            fit.min <- fit
            goodness.of.fit.min <- fit@fit$ics[selection.criteria]
        }
    }
    cat("\n------------------------------------------------------------------------------------------\n")
    cat(paste("Best Model: ", paste(fit.min@formula)[2]))
    cat("\n------------------------------------------------------------------------------------------\n")
    return(fit.min)
}



################################################################################

