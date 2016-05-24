
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
#  .getFormula             Gets the formula.mean and formula.variance from
#                          the object formula.
################################################################################


.getFormula <-
  function(
    formula)
  {
    # Description:
    #   This functions reads formula object and converts it into a list 
    #   containing the separeted mean and variance arguments.
    #   Examples from output:
    #     ~arma(1,1)+garch(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(1,1)
    #     ~aparch(1,1): formula.mean = ~arma(0,0); formula.variance = ~aparch(1,1)
    #     ~arch(1): formula.mean = ~arma(0,0); formula.variance = ~arch(1)
    #     ~arma(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(0,0)
    #     ~ar(1): formula.mean = ~ar(1); formula.variance = ~garch(0,0) 
    #     ~ma(1): formula.mean = ~ma(1); formula.variance = ~garch(0,0) 
    
    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification 
    
    # Return:
    #   A list containing two elements, formula.mean and formula.variance, 
    #   formula.order and the boolean formula.isAPARCH
    
    # FUNCTION: 
    
    # Initial variable declaration
    allLabels = attr(terms(formula), "term.labels")
    formulas.mean.allowed = c("arma")
    formulas.variance.allowed = c("garch","aparch")
    formulaOK <- TRUE
    checkFormulaMean <- ""
    checkFormulaVariance <- ""
    isAPARCH = FALSE
    
    # Error treatment of input parameters 
    if( (length(allLabels) != 1 ) && (length(allLabels) != 2 ) )
      formulaOK <- FALSE
    
    # Formula of type: ~ formula1 + formula2
    else if (length(allLabels) == 2)
    {
      formula.mean = as.formula(paste("~", allLabels[1]))
      formula.var = as.formula(paste("~", allLabels[2]))
      checkFormulaMean = rev(all.names(formula.mean))[1]
      checkFormulaVariance = rev(all.names(formula.var))[1]
      if( !any(formulas.mean.allowed == checkFormulaMean) || 
            !any(formulas.variance.allowed == checkFormulaVariance))   
        formulaOK <- FALSE
    }
    # Formula of type: ~formula1
    else if (length(allLabels) == 1) 
    {
      
      # pure 'garch' or 'aparch'
      if(grepl("arch", attr(terms(formula), "term.labels"))) 
      {
        formula.mean = as.formula("~ arma(0, 0)")
        formula.var = as.formula(paste("~", allLabels[1]))
        checkFormulaVariance = rev(all.names(formula.var))[1]
        if(!any(formulas.variance.allowed == checkFormulaVariance))   
          formulaOK <- FALSE
      }
      else # pure 'ar', 'ma' or 'arma' model.
      {
        formula.mean = as.formula(paste("~", allLabels[1]))  
        formula.var = as.formula("~ garch(0, 0)") 
        checkFormulaMean = rev(all.names(formula.mean))[1]
        if(!any(formulas.mean.allowed == checkFormulaMean))   
          formulaOK <- FALSE
      }    
    }
    
    # Check if we are fitting "aparch" model 
    if(checkFormulaVariance == "aparch")
      isAPARCH = TRUE
    
    # Get model order and check if they were specified correctly
    if(formulaOK == TRUE)
    {
      model.order.mean = 
        as.numeric(strsplit(strsplit(strsplit(as.character(formula.mean), 
                                              "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
      model.order.var = 
        as.numeric(strsplit(strsplit(strsplit(as.character(formula.var), 
                                              "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
      if( (length(model.order.mean) != 2) || (length(model.order.mean) != 2))
        formulaOK <- FALSE 
    }
    
    # Check if model order was specified correctly.
    if(formulaOK == TRUE)
    {
      m = model.order.mean[1]
      n = model.order.mean[2]
      p = model.order.var[1]
      q = model.order.var[2]        
      if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
           any (c(m,n,p,q) < 0) || (p == 0 && q != 0))
        formulaOK <- FALSE            
    }
    
    
    # Stop if formula was not specified correctly
    if(formulaOK == FALSE)
      stop ("Invalid Formula especification. 
            Formula mean must be 'arma' and 
            Formula Variance must be one of: garch or aparch
            For example:
            ARMA(1,1)-GARCH(1,1):  ~arma(1,1)+garch(1,1),
            AR(1)-GARCH(1,1):      ~arma(1,0)+garch(1,1),
            MA(1)-APARCH(1,0):     ~arma(0,1)+aparch(1,0),
            ARMA(1,1):             ~arma(1,1),
            ARCH(2):               ~garch(1,0),
            For more details just type: ?gsFit")
    
    # Return
    list(formula.mean = formula.mean,formula.var = formula.var, 
         formula.order = c(m,n,p,q), isAPARCH = isAPARCH)
  }

################################################################################








