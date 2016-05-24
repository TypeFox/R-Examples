
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



####################################################################################
#  FUNCTION:                        DESCRIPTION:
#
#  .stationarityAparch             Compute the stationarity condition for the 
#                                   GARCH or APARCH model with several cond. 
#                                   distributions
#
#  gsMomentAparch                   Evaluate the exact moments of the type
#                                   E(z - gm|z|)^delta for several conditional
#                                   distributions.
#
#  .normMomentAparch               Exact APARCH moments for the standard normal 
#
#  .stdMomentAparch                Exact APARCH moments for the standard t-Student
#
#  .skstdMomentAparch              Exact APARCH moments for the standard skew
#                                   t-Student defined by Fernandez and Steel (1998).
#                                   Notice that this distribution is standardized 
#                                   with location zero and unit scale, but it is not
#                                   standardized in the sense that it has a zero mean and 
#                                   unit variance. The exact formula for this expression 
#                                   would be very complicated to compute if we reescale 
#                                   the original distribution to (X - mu) / sigma and 
#                                   therefore, we choose to work with the raw distribution
#                                   defined in Fernandez and Steel without this reparametrization
#
#  .gedMomentAparch                Exact APARCH moments for the standard GED distribution

#  .gatMomentAparch                 Exact APARCH moments for the standard GAt distribution
#
#  .stableS1MomentAparch             Exact APARCH moments for the location zero and unit scale
#                                   in S1 parametrization (see Nolan (1999)).
#
#  .stableS1SymmetricMomentGarch    Exact GARCH moments ( E |zt| for the symmetric stable (S1)
#                                   distribution (S1 parametrization)
#
#  .stableS1SymmetricMomentAparch   Exact APARCH moments for the symmetric stable (S1) distribution
#
#  .stableS1MomentPowerGarch        Exact power-GARCH moments for the asymmetric stable (S1) distribution
#
#  .trueAparchMomentsWurtz           Function that evaluate moments with numerical integration.
#                                   This function was used only to test the mathematical 
#                                   formulas implemented here
####################################################################################



.stationarityAparch <-
  function (model = list(), 
            formula,
            cond.dist = c("stableS1", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"))
{
    
    # Description: 
    #   A function that evaluate stationarity conditions 
    #   to guide parameter estimation in GARCH/APARCH models
    #   Make sure we have a garch or aparch with p > 0
    #   The stationarity condition for the APARCH is:
    #   sum ( E[ |z|-gm*z ] ^ delta * alpha + beta ) < 1
    #   and for the GARCH model is
    #   sum ( E[ |z| ] ^ delta * alpha + beta ) < 1
    #   where delta is usually equal to 2, except for the 
    #   stable model where delta = 1.
    #   The expectation  E[ |z| ] is equal to one for the following 
    #   distributions: norm, std, sstd and ged.
    #   For the GAt, gev and stable conditional distribution we need to 
    #   calculate this expectation, even for the GARCH model.  
    
    
    # Arguments:
    #   formula - an object returned by function .getFormula
    #   parm - a list with the model parameters as entries
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gm - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     delta - the exponent value used in the variance equation.
    #     skew - a numeric value listing the distributional
    #        skewness parameter.
    #     shape - a numeric value listing the distributional
    #        shape parameter.
    #   cond.dist - a character string naming the distribution
    #       function.           
    
    # Return Values:
    # sum ( E[ |z|-gm*z ] ^ delta * alpha + beta ) - When it can be computed
    # 1e99 - In case we the expectation E[ |z|-gm*z ] ^ delta could not be evaluated. 
    
    # FUNCTION:
    
    # get parameters
    alpha <- model$alpha
    beta <- model$beta
    delta <- model$delta
    gm <- model$gm
    skew <- model$skew
    shape <- model$shape
    
    # Conditional distribution
    cond.dist = match.arg(cond.dist)
    
   #print(t(model))
    
    
    if(length(alpha) == 0  || length(alpha) != length(gm) || length(delta) != 1)
      stop("Failed to verify conditions:
           if(length(alpha) == 0 || length(delta) != 1)")
    
    # We must have a model with a non zero garch/aparch order
    if(formula$formula.order[3] == 0)
      stop("Invalid model")
    
    # garch model: sum ( E[ |z| ] ^ delta * alpha + beta ) < 1
    if(formula$isAPARCH == FALSE) 
    {
        if (any( c("norm", "std", "sstd", "ged") == cond.dist))
            return(sum(alpha) + sum(beta))  
        
        if(cond.dist == "gev")
            kappa = try(.gevMomentAparch(shape = shape, delta = 2, gm = 0), silent = TRUE)
        
        if(cond.dist == "stableS1")
            kappa = try(.stableS1MomentPowerGarch (shape = shape, skew = skew, 
                                                      delta = 1), silent = TRUE)
        
        if(cond.dist == "gat")
           kappa = try(.gatMomentAparch(shape = shape, skew = skew, 
                                       delta = 2, gm = 0), silent = TRUE)   
        
        if(cond.dist == "skstd")
          kappa = try(.skstdMomentAparch(shape = shape, skew = skew, delta = 2, gm = 0), silent = TRUE)  
        
        if( is.numeric(kappa))
            return(kappa*sum(alpha) + sum(beta)) 
        else 
            return(1e99)      
        
      
    } else {
        # aparch model
        kappa = rep(0,length(alpha))
        
        for( i in 1:length(alpha))
        {
            if(cond.dist == "stableS1")
               kappa[i] = try(.stableS1MomentAparch (shape = shape, skew = skew,
                                                  delta = delta, gm = gm[i]), silent = TRUE)
            if(cond.dist == "gev")
              kappa[i] = try(.gevMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)
            
            if(cond.dist == "gat")
              kappa[i] = try(.gatMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)  
            
            if(cond.dist == "norm")
              kappa[i] = try(.normMomentAparch (delta = delta, gm = gm[i]), silent = TRUE)
            
            if(cond.dist == "std")
              kappa[i] = try(.stdMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)
            
            if(cond.dist == "skstd")
              kappa[i] = try(.skstdMomentAparch(shape = shape, skew = skew, delta = delta, gm = gm[i]), silent = TRUE)
            
            if(cond.dist == "ged")
              kappa[i] = try(.gedMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)     
        }
        if( is.numeric(kappa) )
            return(sum(kappa*alpha) + sum(beta))  
        else 
            return(1e99)
    }   
}



# ------------------------------------------------------------------------------



gsMomentAparch <- function(
    cond.dist = c("stableS1", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"),
    shape = 1.5, 
    skew = 0,
    delta = 1,
    gm = 0)
{
    # conditional distribution
    cond.dist = match.arg(cond.dist)

    
    if(cond.dist == "stableS1")
      kappa = .stableS1MomentAparch (shape = shape, skew = skew,
                                           delta = delta, gm = gm)
    if(cond.dist == "gev")
      kappa = .gevMomentAparch(shape = shape, delta = delta, gm = gm)
    
    if(cond.dist == "gat")
      kappa = .gatMomentAparch(shape = shape, delta = delta, skew = skew, gm = gm)  
    
    if(cond.dist == "norm")
      kappa = .normMomentAparch (delta = delta, gm = gm)
    
    if(cond.dist == "std")
      kappa = .stdMomentAparch(shape = shape, delta = delta, gm = gm)
    
    if(cond.dist == "sstd")
      kappa = .sstdMomentAparch(shape = shape, skew = skew, delta = delta, gm = gm)
    
    if(cond.dist == "skstd")
      kappa = .skstdMomentAparch(shape = shape, skew = skew, delta = delta, gm = gm)
    
    if(cond.dist == "ged")
      kappa = .gedMomentAparch(shape = shape, delta = delta, gm = gm)
    
    # Return
    kappa

}



# ------------------------------------------------------------------------------



.normMomentAparch <- function(delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard normal distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: A long memory property of stock market returns and a 
  #   new model. Ding, Granger and Engle (1993), Appendix B. 
  
  # Error treatment of input parameters
  if( ( abs ( gm ) >= 1 ) || ( delta <= 0 ) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following condition cannot be true
         ( abs ( gm ) >= 1 ) || ( delta <= 0 )")
  
  result = 1 / sqrt ( 2 * pi ) * ( (1 + gm ) ^ delta + ( 1 - gm ) ^ delta ) *
    2 ^ ( ( delta - 1 ) / 2 ) * gamma ( ( delta + 1 ) / 2 )
  
  # Return
  result
}



# ------------------------------------------------------------------------------


.stdMomentAparch <- function(shape = 3, delta = 1.2, gm = 0)
{
    # Description:
    #   Returns the following Expectation for a standard t-Student distribution
    #   E[ (|z|-gm*z)^delta.
    #   Reference: Mittnik -  (2000) - Conditional Density and Value-at-Risk 
    #   Prediction of Asian Currency Exchange Rates
    #   Note that the formula provided by Mittnik is valid for the original t-Student 
    #   distribution (which is nonStandardized).
    #   Since our objective is to calculate it for the standart t-Student distribution
    #   defined in Wurtz et al. (2006) - eq (13) we need to multiply the Mittnik
    #   formula by (sqrt((shape-2)/shape))^delta.
    #   The main drawback is that the standardized t-Student distribution 
    #   is valid only for shape > 2 ( in which case the variance is finite)
    # Error treatment of input parameters
    if( (abs(gm) >= 1) || (delta <= 0) || (shape <= 2) || (delta >= shape))
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           (abs(gm) >= 1) || (delta <= 0) || (shape <= 0)|| (delta >= shape)")
    
    mittnikFormula = shape^(delta/2)*1/2/sqrt(pi)*( (1+gm)^delta+(1-gm)^delta )*
      gamma((delta+1)/2)*(gamma(shape/2))^(-1)*gamma((shape-delta)/2)
    factorToMultiply = (sqrt((shape-2)/shape))^delta
    
    # Return
    mittnikFormula*factorToMultiply
}


# ------------------------------------------------------------------------------

.skstdMomentAparch <- function(shape = 3, skew = 1, delta = 1.2, gm = 0)
{
    # Description:
    #   Returns the following Expectation for a standard skew t-Student distribution
    #   E[ (|z|-gm*z)^delta.
    #   Reference: fGarch: Wurtz and the package Skewt uses the one defined 
    #   in Fernandez, C. and Steel, M. F. J. (1998). On Bayesian modeling of fat 
    #   tails and skewness, J. Am. Statist. Assoc. 93, 359-371.
    #   Another reference: http://www.timberlake.co.uk/slaurent/G@RCH/Book63.html
     
    # Note: This is not the APARCH moments for the 'dsstd' function from fGarch package. 
    # This is the skew t-Student defined by Fernandez and Steel without the 
    # reparametrization introduced by Wurtz et al. (2006).
    # but theThe point here is that the SPLUS FinMetrics
    # also uses this distribution without reparametrization. 
    # Indeed, Mittnik et al. (2000) also uses his t3-distribution (GAt) without such
    # reparametrization. These distributions are in fact standardized, but the 
    # point is that the assymetry parameter plays a crucial role in the calculation
    # of the APARCH moments.   
    
    if( (abs(gm) >= 1) || (delta <= 0) || (shape <= 2) || 
        (delta >= shape) || (skew <= 0)) 
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           (abs(gm) >= 1) || (delta <= 0) || (shape <= 2) || 
           (delta >= shape) || (skew <= 0)")
    a = skew^(-1-delta)*(1+gm)^delta + skew^(1+delta)*(1-gm)^delta
    b = gamma((delta+1)/2)*gamma((shape-delta)/2)*(shape-2)^((1+delta)/2)
    c = (skew + 1/skew)*sqrt((shape-2)*pi)*gamma(shape/2)
    
    # Return
    a*b/c  
}




# ------------------------------------------------------------------------------




.gatMomentAparch <- function(shape = c(3,1), skew = 1, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard GAt distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: Mittnik -  (2000) - Conditional Density and Value-at-Risk 
  #   Prediction of Asian Currency Exchange Rates.
  
  if( (abs(gm) >= 1) || (delta <= 0) || (shape[1] <= 0) || (shape[2] <= 0) ||
        (delta >= (shape[1])*(shape[2]) ) || (skew <= 0)) 
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true
         (abs(gm) >= 1) || (delta <= 0) || (shape[1] <= 0) || (shape[2] <= 0)
         (delta >= shape[1]*shape[2]) || (skew <= 0)")
  
    nu = shape[1]
    d = shape[2]
    theta = skew
  
    a = ( 1 + gm ) ^ delta * theta ^ ( - delta - 1 ) + ( 1 - gm ) ^ delta * theta ^ ( delta + 1 )
    b = nu ^ ( delta / d ) * beta ( ( delta + 1 ) / d , nu - delta / d )
    c = ( 1 / theta + theta ) * beta( 1 / d , nu )
  
    # Return
    a * b / c  
}




# ------------------------------------------------------------------------------


.gedMomentAparch <- function(shape = 3, delta = 1.2, gm = 0)
{
    # Description:
    #   Returns the following Expectation for a standard skew t-Student distribution
    #   E[ (|z|-gm*z)^delta.
    #   Reference: http://www.timberlake.co.uk/slaurent/G@RCH/Book63.html
    
    if( (abs(gm) >= 1) || (delta <= 0)) 
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
             The following conditions cannot be true
             ")
    lambda = sqrt(gamma(1/shape)*2^(-2/shape)/gamma(3/shape))
    ((1+gm)^delta + (1-gm)^delta)*2^((delta-shape)/shape)*
    gamma((delta+1)/shape)*lambda^delta/gamma(1/shape)
}



# ------------------------------------------------------------------------------



.gevMomentAparch <- function(shape = 0.3, delta = 1.2, gm = 0)
{
  
    # FIND A BETTER WAY TO CALCULATE IT INSTEAD OF USING THE INTEGRATION FUNCTION
  
  
    # Description:
    #   Returns the following Expectation for a standard (location zero and scale 1) 
    #   GEV distribution 
    #   E[ (|z|-gm*z)^delta.
    #   Reference: GEVStableGarch papper
    
    if( ( abs( gm ) >= 1 ) || ( delta <= 0 )) 
    # if( ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= -0.5 ) ) 
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           ( abs( gm ) >= 1 ) || ( delta <= 0 )")
    
    xi = shape
    
    e = function(x, gm, delta) {
      (abs(x)-gm*x)^delta * dgev(x, xi = xi)
    }
    
    # Compute Moment by Integration
    I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.25,
                  gm = gm, delta = delta)
    
    # Return
    as.numeric(I[1])
}

# ------------------------------------------------------------------------------



.sstdMomentAparch <- function(shape = 4, skew = 1, delta = 1.2, gm = 0)
{
  
    # FIND A BETTER WAY TO CALCULATE IT INSTEAD OF USING THE INTEGRATION FUNCTION
    
    
    # Description:
    #   Returns the following Expectation for the sstd distribution implemented 
    #   inside package fGarch
    #   sstd distribution
    #   E[ (|z|-gm*z)^delta.
    #   Reference: fGarch package
    
    if( ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= 2 ) || ( skew <= 0 )) 
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
             The following conditions cannot be true
             ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= 2 ) || ( skew <= 0 )")
    
    nu = shape
    xi = skew
    
    e = function(x, gm, delta) {
      (abs(x)-gm*x)^delta * dsstd(x, nu = nu, xi = xi)
    }
    
    # Compute Moment by Integration
    I = integrate(e, lower = -Inf, upper = Inf, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.25,
                  gm = gm, delta = delta)
    
    # Return
    as.numeric(I[1])
}



# ------------------------------------------------------------------------------



.stableS1MomentAparch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
  # Description: 
  #   Returns the following Expectation for a standard stable distribution in S1 parametrization
  #   E[ (|z|-gm*z)^delta.
  #   This formula uses S(alpha,skew,1,0;pm) where pm = 1.
  #   Reference: The GEVStableGarch papper on JSS. 
  
  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
      (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
         (delta <= 1) || (delta >= shape)
         Note: This function do not accept the normal case, i.e, alpha = 2")
  
  # Calculate the expression
  sigma.til <- (1 + (skew*tan(shape*pi/2))^2)^(1/2/shape)
  k.shape <- shape - 2 # since skew > 1
  skew.til <- 2/pi/(shape-2)*atan(skew*tan((shape-2)*pi/2))
  g1 <- gamma((1/2 + skew.til*k.shape/2/shape)*(-delta))
  g2 <- gamma(1/2 - skew.til*k.shape/2/shape + (1/2 + skew.til*k.shape/2/shape)*(delta+1))
  g3 <- gamma((1/2 - skew.til*k.shape/2/shape)*(-delta))
  g4 <- gamma(1/2 + skew.til*k.shape/2/shape + (1/2 - skew.til*k.shape/2/shape)*(delta+1))
  kappa <- 1/shape/sigma.til*sigma.til^(delta+1)*gamma(delta+1)*gamma(-delta/shape)*
    (1/g1/g2*(1-gm)^delta + 1/g3/g4*(1+gm)^delta)
  
  # Return 
  kappa
}









# ------------------------------------------------------------------------------





.stableS0MomentAparch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
  # Description: 
  #   Returns the following Expectation for standard stable distribution in S0 parametrization
  #   E[ (|z|-gm*z)^delta.
  #   This formula uses S(alpha,skew,1,0;pm) where pm = 0.
  #   Reference: The GEVStableGarch papper on JSS. 
  
  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
      (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
         (delta <= 1) || (delta >= shape)
         Note: This function do not accept the normal case, i.e, alpha = 2")
  
  # Compute Moment by Integration
    e = function(x, gm, delta) {
      (abs(x)-gm*x)^delta * stable::dstable.quick(x = x, alpha = shape,
                      beta = skew, param = 0)
    }
	
    I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.25,
                  gm = gm, delta = delta)
    
    # Return
    as.numeric(I[1])
}


.stableS2MomentAparch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
  # Description: 
  #   Returns the following Expectation for standard stable distribution in S2 parametrization
  #   E[ (|z|-gm*z)^delta.
  #   This formula uses S(alpha,skew,1,0;pm) where pm = 2.
  #   Reference: The GEVStableGarch papper on JSS. 
  
  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
      (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
         (delta <= 1) || (delta >= shape)
         Note: This function do not accept the normal case, i.e, alpha = 2")
  
  # Compute Moment by Integration
    e = function(x, gm, delta) {
      (abs(x)-gm*x)^delta * stable::dstable.quick(x = x, alpha = shape,
                      beta = skew, param = 2)
    }
	
    I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.25,
                  gm = gm, delta = delta)
    
    # Return
    as.numeric(I[1])
}






.gevMomentAparch <- function(shape = 0.3, delta = 1.2, gm = 0)
{
  
    # FIND A BETTER WAY TO CALCULATE IT INSTEAD OF USING THE INTEGRATION FUNCTION
  
  
    # Description:
    #   Returns the following Expectation for a standard (location zero and scale 1) 
    #   GEV distribution 
    #   E[ (|z|-gm*z)^delta.
    #   Reference: GEVStableGarch papper
    
    if( ( abs( gm ) >= 1 ) || ( delta <= 0 )) 
    # if( ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= -0.5 ) ) 
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           ( abs( gm ) >= 1 ) || ( delta <= 0 )")
    
    xi = shape
    
    e = function(x, gm, delta) {
      (abs(x)-gm*x)^delta * dgev(x, xi = xi)
    }
    
    # Compute Moment by Integration
    I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.25,
                  gm = gm, delta = delta)
    
    # Return
    as.numeric(I[1])
}

# ------------------------------------------------------------------------------










# ------------------------------------------------------------------------------
.stableS1SymmetricMomentGarch <- function(shape = 1.5)
{
  # See Mittinik et al. (1995) for the definition of 
  # the stable garch model with conditional symmetric stable distribution
  # This formula uses S(alpha,skew,1,0;pm) where pm = 1. 
  if( (shape <= 1) || (shape > 2) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 0) || (shape > 2)")
  
  # Return
  gamma(1 - 1/shape)*2/pi
}



# ------------------------------------------------------------------------------


.stableS1SymmetricMomentAparch <- function(shape = 1.5, delta = 1.2, gm = 0)
{
  # See Diongue - 2008 (An investigation of the stable-Paretian Asymmetric Power GARCH model)
  # This formula uses S(alpha,0,1,0;pm) where pm = 1. 
  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true.
            (shape <= 1) || (shape > 2) || (delta >= shape))")
  
  
  # k(delta): At this point we have delta > 1
  if(delta == 1)
    k = pi/2
  else
    k = gamma(1-delta)*cos(pi*delta/2)
  result = k^(-1)*gamma(1-delta/shape)*1/2*((1+gm)^delta+(1-gm)^delta)
  
  # Return
  result
}







# ------------------------------------------------------------------------------
.stableS1MomentPowerGarch <- function(shape = 1.5, skew = 0.5, delta = 1.2)
{
    # See Mittnik et al. (2002) - Stationarity of stable power-GARCH processes
    # This formula uses S(alpha,skew,1,0;pm) where pm = 1 (tests reveald) 
    # Error treatment of input parameters
    if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || 
        (delta >= shape))
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true.
            (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || 
            (delta >= shape)")
    
    
    # k(delta): At this point we have delta > 1
    if(delta == 1)
        k = pi/2
    else
        k = gamma(1-delta)*cos(pi*delta/2)
    # tal(shape,skew)
    tau = skew*tan(shape*pi/2)
    # lambda(shape,skew,delta): eq (10) Mittnik et al. (2002)
    lambda = k^(-1)*gamma(1-delta/shape)*(1 + tau^2)^(delta/2/shape)*cos(delta/shape*atan(tau))
    lambda
}



# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
.trueAparchMomentsWurtz <- 
  function(fun = "norm", gm = 0, delta = 1, lower = -Inf, upper = Inf, ...)
  { 
    # This function was originally from package fGarch. We 
    # changed it to return only the aparch moment expression.
    # Description:
    #   Computes APARCH moments using Integration
    
    # Arguments:
    #   fun - name of density functions of APARCH innovations
    #   alpha, gm - numeric value or vector of APARCH coefficients,
    #       must be of same length  
    #   beta - numeric value or vector of APARCH coefficients
    #   delta - numeric value of APARCH exponent
    
    # Note:
    #   fun is one of: norm, snorn, std, sstd, ged, sged, snig
    
    # FUNCTION:
    
    # Match Density Function:
    fun = match.fun(fun)
    
    # Persisgtence Function: E(|z|-gm z)^delta
    e = function(x, gm, delta, ...) {
      (abs(x)-gm*x)^delta * fun(x, ...)
    }
    
    # Compute Moment by Integration
    I = integrate(e, lower = lower, upper = upper, subdivisions = 1000,
                    rel.tol = .Machine$double.eps^0.25,
                    gm = gm, delta = delta, ...)
    
    # Return Value:
    I
  }

################################################################################