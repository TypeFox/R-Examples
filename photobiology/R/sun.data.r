#' @title Solar spectral irradiance (simulated)
#'
#' @description A dataset containing the wavelengths at a 1 nm interval and the
#'   corresponding spectral (energy) irradiance and spectral photon irradiance.
#'   Values simulated for 22 June 2010, near midday, at Helsinki, under partly
#'   cloudy conditions. The variables are as follows:
#'
#' @details \itemize{ \item w.length (nm), range 293 to 800 nm. \item s.e.irrad
#' (W m-2 nm-1) \item s.q.irrad (mol m-2 nm-1) }
#'
#'
#' @author Anders K. Lindfors (data)
#' @references Lindfors, A.; Heikkilä, A.; Kaurola, J.; Koskela, T. & Lakkala,
#' K. (2009) Reconstruction of Solar Spectral Surface UV Irradiances Using
#' Radiative Transfer Simulations. Photochemistry and Photobiology, 85:
#' 1233–1239
#'
#' @docType data
#' @keywords datasets
#' @format A \code{source_spct} object with 508 rows and 3 variables
#' @name sun.spct
NULL

#' @title Solar spectral irradiance (simulated)
#'
#' @description A dataset containing the wavelengths at a 1 nm interval and the
#'   corresponding spectral (energy) irradiance and spectral photon irradiance.
#'   Values simulated for 22 June 2010, near midday, at Helsinki, under partly
#'   cloudy conditions. The variables are as follows:
#'
#' @details \itemize{ \item w.length (nm), range 293 to 800 nm. \item s.e.irrad
#' (W m-2 nm-1) \item s.q.irrad (mol m-2 nm-1) }
#'
#' @author Anders K. Lindfors (data)
#'
#' @references Lindfors, A.; Heikkilä, A.; Kaurola, J.; Koskela, T. & Lakkala,
#' K. (2009) Reconstruction of Solar Spectral Surface UV Irradiances Using
#' Radiative Transfer Simulations. Photochemistry and Photobiology, 85:
#' 1233–1239
#'
#' @docType data
#' @keywords datasets
#' @format A \code{data.frame} object with 508 rows and 3 variables
#' @name sun.data
NULL

#' Daily solar spectral irradiance (simulated)
#'
#' A dataset containing the wavelengths at a 1 nm interval and the corresponding
#' spectral (energy) irradiance. Values simulated for 2 June 2012, at Helsinki,
#' under clear sky conditions. The variables are as follows:
#'
#' \itemize{ \item w.length (nm), range 290 to 800 nm. \item s.e.irrad (J d-1
#' m-2 nm-1) \item s.q.irrad (mol d-1 m-2 nm-1) }
#'
#' @author Anders K. Lindfors (data)
#' @references Lindfors, A.; Heikkilä, A.; Kaurola, J.; Koskela, T. & Lakkala,
#' K. (2009) Reconstruction of Solar Spectral Surface UV Irradiances Using
#' Radiative Transfer Simulations. Photochemistry and Photobiology, 85:
#' 1233–1239
#'
#' @note The simualtions are based on libRadTran using hourly mean global
#' radiation measurements to estimate cloud cover. The simulations were for
#' each hour and the results integrated for the whole day.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{source_spct} object with 511 rows and 3 variables
#' @name sun.daily.spct
NULL

#' Daily solar spectral irradiance (simulated)
#'
#' A dataset containing the wavelengths at a 1 nm interval and the corresponding
#' spectral (energy) irradiance. Values simulated for 2 June 2012, at Helsinki,
#' under clear sky conditions. The variables are as follows:
#'
#' \itemize{ \item w.length (nm), range 290 to 800 nm. \item s.e.irrad (J d-1
#' m-2 nm-1) \item s.q.irrad (mol d-1 m-2 nm-1) }
#'
#' @author Anders K. Lindfors (data)
#' @references Lindfors, A.; Heikkilä, A.; Kaurola, J.; Koskela, T. & Lakkala,
#' K. (2009) Reconstruction of Solar Spectral Surface UV Irradiances Using
#' Radiative Transfer Simulations. Photochemistry and Photobiology, 85:
#' 1233–1239
#'
#' @docType data
#' @keywords datasets
#' @format A \code{data.frame} object with 511 rows and 3 variables
#' @name sun.daily.data
NULL
