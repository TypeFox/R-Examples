AllModuleStatistics3D <- function (modules, ndig = 4) 
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

    n_elements = NumberOfElements3D(modules)
	min_value = MinValue3D(modules)
    max_value = MaxValue3D(modules)
    range_value = Range3D(modules)
    module_sum = ModuleSum3D(modules)
    m_arithmetic = ArithmeticMean3D(modules)
    s_error = StandardError3D(modules)
    s_d_module = ModuleStandardDeviation3D(modules)
    v_module = ModuleVariance3D(modules)
    s_d_module_p = ModulePopulationStandardDeviation3D(modules)
    v_module_p = ModulePopulationVariance3D(modules)
	cs = SkewnessModuleCoefficient3D(modules)
    ca = KurtosisModuleCoefficient3D(modules)
	
    print("  ---------------------------  ")
    print("  LINEAR STATISTICS - MODULES  ")
    print("  ---------------------------  ")	
    print(paste("  NUMBER OF ELEMENTS =", format(round(n_elements, 0), nsmall = 0)))
    print(paste("  MIN VALUE =", format(round(min_value, ndig), nsmall = ndig)))
    print(paste("  MAX VALUE =", format(round(max_value, ndig), nsmall = ndig)))
    print(paste("  RANGE =", format(round(range_value, ndig), nsmall = ndig)))
    print(paste("  ARITHMETIC MEAN =", format(round(m_arithmetic, 4), nsmall = ndig)))
    print(paste("  MEAN STANDARD ERROR =", format(round(s_error, ndig), nsmall = ndig)))
    print(paste("  STANDARD DEVIATION =", format(round(s_d_module, ndig), nsmall = ndig)))
    print(paste("  VARIANCE = ", format(round(v_module, ndig), nsmall = ndig)))
    print(paste("  POPULATION STANDARD DEVIATION =", format(round(s_d_module_p, ndig), nsmall = ndig)))   
    print(paste("  POPULATION VARIANCE =", format(round(v_module_p, ndig), nsmall = ndig)))
    print(paste("  SKEWNESS COEFFICIENT =", format(round(cs, ndig), nsmall = ndig)))
    print(paste("  KURTOSIS COEFFICIENT =", format(round(ca, ndig), nsmall = ndig)))
}