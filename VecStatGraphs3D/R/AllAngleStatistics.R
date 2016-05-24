AllAngleStatistics <- function (incr, ndig = 4) 
{
	
  # Calculates the basic spherical statistics from the 
  # X, Y and Z coordinate increments
  #
  # Args:
  #   coord : values of the coordinate increments (matrix)
  #   ndig  : decimal places
  #
  # Returns:
  #   Sample size
  #   Mean colatitude
  #   Mean longitude
  #   Mean module
  #   Concentration parameter
  #   Spherical standard error
  #

	n_elements = dim(incr)[1]
	
	vm_direction <- MeanDirection3D(incr)
    vm_module <- MeanModule3D(incr)
	
	unit_incr <- RealModToUnitMod3D(incr)
    um_direction <- MeanDirection3D(unit_incr)
    um_module <- MeanModule3D(unit_incr)
    conc_parameter <- ConcentrationParameter3D(unit_incr)
	sphericalErr <- SphericalStandardError3D(unit_incr)
	
    print("  -------------------------------  ")
    print("  SPHERICAL STATISTICS - ANGLES    ")
    print("  -------------------------------  ")
	print(paste("  NUMBER OF ELEMENTS =", format(round(n_elements, 0), nsmall = 0)))
	print("                                   ")
	print("  Statistics for real (non-unit) vectors  ")
    print("  --------------------------------------  ")
    print(paste("  COLATITUDE  =", format(round(vm_direction[1], ndig), nsmall = ndig)))
 	print(paste("  LONGITUDE   =", format(round(vm_direction[2], ndig), nsmall = ndig)))
    print(paste("  MEAN MODULE =", format(round(vm_module, ndig), nsmall = ndig)))
	print("                                 ")
	print("  Statistics for unit vectors    ")
    print("  -----------------------------  ")
	print(paste("  COLATITUDE  =", format(round(um_direction[1], ndig), nsmall = ndig)))
 	print(paste("  LONGITUDE   =", format(round(um_direction[2], ndig), nsmall = ndig)))
    print(paste("  MEAN MODULE =", format(round(um_module, ndig), nsmall = ndig)))
    print(paste("  CONCENTRATION PARAMETER  =", format(round(conc_parameter, ndig), nsmall = ndig)))
    print(paste("  SPHERICAL STANDARD ERROR =", format(round(sphericalErr, ndig), nsmall = ndig)))
}
