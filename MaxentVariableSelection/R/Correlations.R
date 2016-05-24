Correlations <-
    function(important.variable,variablenames,backgroundsites,correlationthreshold){
                                        # variables that are correlated with a correlation coefficient higher than the 'correlationthreshold' will be  excluded
        
                                        # The first variablename is expected to be the variable with highest contribution to the maxent model
        
        
        backgroundvalues <- read.csv(backgroundsites,sep=",",header=TRUE)
                                        # select the specified variables from the backgroundvalues
        backgroundvalues <- backgroundvalues[,(colnames(backgroundvalues) %in% variablenames)]
        
                                        # calculate correlation coefficients
        correlation.coefficients <- cor(backgroundvalues)
        
                                        # Extract the correlation coefficients to the most important variable
        correlation.coefficients <- correlation.coefficients[,colnames(correlation.coefficients)==important.variable]
        
        uncorrelatedvariables <-  c(important.variable,names(correlation.coefficients)[abs(correlation.coefficients)<correlationthreshold])
        return(list(correlation.coefficients,uncorrelatedvariables))
    }
