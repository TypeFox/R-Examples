#' @details
#' This package provides constructors for objects of class \code{waveband}
#' from package 'photobiology'. These contructors are based on standard
#' definitions and frequently used non-standardized definitions. When
#' different definitions are in common use for a given named waveband the
#' contructors accept an argument to chose among them.
#'
#' By necesity we cover only a subset of all definitions in use. These should
#' be thought as convenience functions, as waveband objects according to any
#' arbitrary definition can be constructed with the functions provided by
#' package \code{\link[photobiology]{photobiology-package}}
#'
#' @references
#' Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A. R., Robson, T. M.,
#' Rosenqvist, E. (Eds.). (2012). Beyond the Visible: A handbook of best
#' practice in plant UV photobiology (1st ed., p. xxx + 174). Helsinki:
#' University of Helsinki, Department of Biosciences, Division of Plant Biology.
#' ISBN 978-952-10-8363-1 (PDF), 978-952-10-8362-4 (paperback). Open access PDF
#' download available at http://hdl.handle.net/10138/37558
#'
#' Caldwell, M. M. (1971) Solar UV irradiation and the growth and development
#' of higher plants. In Giese, A. C. (Ed.) Photophysiology, Academic Press,
#' 1971, 6, 131-177
#'
#' Diffey, B. L. 1991. Solar ultraviolet radiation effects on biological systems.
#' Review in Physics in Medicine and Biology 36 (3): 299-328.
#'
#' Green, A. E. S., Miller, J. H. (1975) Measures of biologically active
#' radiation in the 280-340 nm region. Impacts of climate change on the
#' environment. CIAP Monograph, 5, Part 1, Chapter 2.2.4
#'
#' Green, A. E. S.; Sawada, T. & Shettle, E. P. (1974) The middle
#' ultraviolet reaching the ground Photochemistry and Photobiology, 1974, 19,
#' 251-259
#'
#' Ibdah, M., Krins, A., Seidlitz, H. K., Heller, W., Strack, D. & Vogt, T.
#' (2002) Spectral dependence of flavonol and betacyanin accumulation in
#' Mesembryanthemum crystallinum under enhanced ultraviolet radiation. Plant,
#' Cell & Environment, 25, 1145-1154
#'
#' INTERNATIONAL COMMISSION ON NON-IONIZING RADIATION PROTECTION (2004) ICNIRP
#' GUIDELINES ON LIMITS OF EXPOSURE TO ULTRAVIOLET RADIATION OF WAVELENGTHS
#' BETWEEN 180 nm AND 400 nm (INCOHERENT OPTICAL RADIATION). HEALTH PHYSICS
#' 87(2):171-186. \url{http://www.icnirp.org/cms/upload/publications/ICNIRPUV2004.pdf}
#'
#' ISO (2007) Space environment (natural and artificial) - Process for
#' determining solar irradiances. ISO Standard 21348. ISO, Geneva.
#'
#' Quaite, F. E., Sutherland, B. M., Sutherland, J. C. Action spectrum for DNA
#' damage in alfalfa lowers predicted impact of ozone depletion. Nature, 1992,
#' 358, 576–578
#'
#' Micheletti, M. I.; Piacentini, R. D. & Madronich, S. (2003) Sensitivity
#' of Biologically Active UV Radiation to Stratospheric Ozone Changes: Effects
#' of Action Spectrum Shape
#' and Wavelength Range Photochemistry and Photobiology, 78, 456-461
#'
#' Musil, C. F. (1995) Differential effects of elevated ultraviolet-B radiation on the
#' photochemical and reproductive performances of dicotyledonous and
#' monocotyledonous arid-environment ephemerals Plant, Cell and Environment,
#' 18, 844-854
#'
#' Murakami, K., Aiga I. (1994) Red/Far-red photon flux ratio used as
#' an index number for morphological control of plant growth under
#' artificial lighting conditions. Proc. Int. Symp. Artificial Lighting,
#' Acta Horticulturae, 418, ISHS 1997.
#'
#' Sellaro, R., Crepy, M., Trupkin, S. A., Karayekov, E., Buchovsky, A. S.,
#' Rossi, C., & Casal, J. J. (2010). Cryptochrome as a sensor of the blue/green
#' ratio of natural radiation in Arabidopsis. Plant physiology, 154(1), 401-409.
#' doi:10.1104/pp.110.160820
#'
#' Setlow, R. B. (1974) The Wavelengths in Sunlight Effective in Producing Skin
#' Cancer: A Theoretical Analysis. Proceedings of the National Academy of
#' Sciences, 71, 3363-3366
#'
#' Smith, H. (1982) Light quality, photoperception and plant strategy. Annual
#' Review of Plant Physiology, 33:481-518.
#'
#' Webb, A. R.; Slaper, H.; Koepke, P. & Schmalwieser, A. W. Know your standard:
#' clarifying the CIE erythema action spectrum. Photochemistry and photobiology,
#' 2011, 87, 483-486
#'
#' @import photobiology
#' @examples
#'
#' q_irrad(sun.spct, PAR())  # PAR photon irradiance
#' q_irrad(sun.spct, Blue("ISO")) # blue photon irradiance, ISO definition
#' q_irrad(sun.spct, Blue("Sellaro")) # blue photon irradiance, Sellaro et al.'s definition
#' e_irrad(sun.spct, VIS()) # VIS irradiance, ISO definition
#' q_irrad(sun.spct, VIS()) # VIS photon, ISO definition
"_PACKAGE"
