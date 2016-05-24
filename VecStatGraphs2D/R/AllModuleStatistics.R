AllModuleStatistics <- function (modules, ndig = 4) 
{

  # Calculates the basic linear statistics of a sample of linear magnitudes
  #
  # Args:
  #   modules : modules or another linear magnitude
  #   ndig : decimal places
  #
  # Returns:
  #   Sample size
  #   Min, max and range values
  #   Arithmetic mean
  #   Standard error of the mean 
  #   Sample and population standard deviation
  #   Sample and population variance
  #   Skewness coefficient
  #   Kurtosis coefficient
  #

    n_elements = length(modules) 
    max_value = max(modules)
    min_value = min(modules)
    range_value = abs(max_value - min_value)
    m_arithmetic = sum(modules) / length(modules)			
    s_error = StandardError(modules)
    s_d_module = ModuleStandardDeviation(modules)
    v_module = ModuleVariance(modules)
    s_d_module_p = ModulePopulationStandardDeviation(modules)
    v_module_p = ModulePopulationVariance(modules)
    cs = SkewnessModuleCoefficient(modules)
    ca = KurtosisModuleCoefficient(modules)
	
    print("  ---------------------------  ")
    print("  LINEAR STATISTICS - MODULES  ")
    print("  ---------------------------  ")
    print(paste("  NUMBER OF ELEMENTS = ", format(round(n_elements, 0), nsmall = 0)))
    print(paste("  MIN VALUE = ", format(round(min_value, ndig), nsmall = ndig)))
    print(paste("  MAX VALUE = ", format(round(max_value, ndig), nsmall = ndig)))
    print(paste("  RANGE = ", format(round(range_value, ndig), nsmall = ndig)))
    print(paste("  ARITHMETIC MEAN = ", format(round(m_arithmetic, 4), nsmall = ndig)))
    print(paste("  MEAN STANDARD ERROR = ", format(round(s_error, ndig), nsmall = ndig)))
    print(paste("  STANDARD DEVIATION = ", format(round(s_d_module, ndig), nsmall = ndig)))
    print(paste("  VARIANCE = ", format(round(v_module, ndig), nsmall = ndig)))
    print(paste("  POPULATION STANDARD DEVIATION = ", format(round(s_d_module_p, ndig), nsmall = ndig)))   
    print(paste("  POPULATION VARIANCE = ", format(round(v_module_p, ndig), nsmall = ndig)))
    print(paste("  SKEWNESS COEFFICIENT = ", format(round(cs, ndig), nsmall = ndig)))
    print(paste("  KURTOSIS COEFFICIENT = ", format(round(ca, ndig), nsmall = ndig)))
    
}
