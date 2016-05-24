.validate <- function (object) {
            ncol <- ncol (object@data$spc)
				
            if (is.null (ncol))
              ncol <- 0
            
			  	if (length (object@wavelength) != ncol)
              return ("Length of wavelength vector differs from number of data points per spectrum.")

            TRUE
          }
