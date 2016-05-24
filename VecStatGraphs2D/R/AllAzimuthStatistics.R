AllAzimuthStatistics <- function (azimuths, ndig = 4) 

{

  # Calculates the basic circular statistics from a sample of azimuths
  #
  # Args:
  #   azimuths : azimuths (degrees)
  #   ndig : decimal places
  #
  # Returns:
  #   Sample size
  #   Mean azimuth
  #   Mean module
  #   Circular standard deviation
  #   Circular variance
  #   Circular dispersal
  #   Von Mises parameter
  #   Skewness coefficient
  #   Kurtosis coefficient
  #

    n_elements = length(azimuths)
    m_azimuth = MeanAzimuth(azimuths)
    m_module = MeanModule(azimuths)
    c_variance = CircularVariance(azimuths)
    s_deviation = CircularStandardDeviation(azimuths)
    vm_parameter = VonMisesParameter(azimuths)
    c_dispersal = CircularDispersal(azimuths)
    s_azimuth = SkewnessAzimuthCoefficient(azimuths)
    k_azimuth = KurtosisAzimuthCoefficient(azimuths)
    
    print("  ------------------------------  ")
    print("  CIRCULAR STATISTICS - AZIMUTHS  ")
    print("  ------------------------------  ")
    print(paste("  NUMBER OF ELEMENTS =", format(round(n_elements, 0), nsmall = 0)))
    print(paste("  MEAN AZIMUTH =", format(round(m_azimuth, ndig), nsmall = ndig)))
    print(paste("  MEAN MODULE =", format(round(m_module, ndig), nsmall = ndig)))
    print(paste("  CIRCULAR STANDARD DEVIATION =", format(round(s_deviation, ndig), nsmall = ndig)))
    print(paste("  CIRCULAR VARIANCE =", format(round(c_variance, ndig), nsmall = ndig)))
    print(paste("  CIRCULAR DISPERSAL =", format(round(c_dispersal, ndig), nsmall = ndig)))
    print(paste("  VON MISES PARAMETER =", format(round(vm_parameter, ndig), nsmall = ndig)))
    print(paste("  SKEWNESS COEFFICIENT =", format(round(s_azimuth, ndig), nsmall = ndig)))
    print(paste("  KURTOSIS COEFFICIENT =", format(round(k_azimuth, ndig), nsmall = ndig)))

}
