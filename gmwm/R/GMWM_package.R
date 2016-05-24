# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' Generalized Method of Wavelet Moments (GMWM) Package
#'
#' Generalized Method of Wavelet Moments (GMWM) is an estimation technique for the parameters of time series models. It uses the wavelet variance in a moment matching approach that makes it particularly suitable for the estimation of certain state-space models. Furthermore, there exists a robust implementation of GMWM, which allows the robust estimation of some state-space models and ARIMA models. Lastly, the package provides the ability to quickly generate time series data, perform different wavelet decompositions, and visualizations. 
#' 
#' @details 
#' \tabular{ll}{
#' Package: \tab GMWM\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.1\cr
#' Date: \tab 2015-12-28\cr
#' License: \tab CC BY-NC-SA 4.0\cr
#' }
#' 
#' @section Financial Support:
#' 
#' Primary support for this project was made possible by the Department of Statistics at the University of Illinois at Urbana-Champaign (UIUC). Furthermore, the contributions of undergraduate researcher Wenchao Yang for improving the visualization systems was supported, in part, by a grant from the Office of Undergraduate Research at UIUC. The project was mentored by Prof. Stephane Guerrier. 
#' 
#' @section Acknowledgements:
#' 
#' We would like to thank the following individuals for their contributions and advice in the development of the GMWM methodology:
#' 
#' \itemize{
#'  \item Prof. Maria-Pia Victoria-Feser 
#'  \item Dr. Jan Skaloud
#'  \item Dr. Yannick Stebler
#' }
#' 
#' Furthermore, we are also greatful to Dr. Jan Skaloud and Philipp Clausen of Geodetic Engineering Laboratory (TOPO), Swiss Federal Institute of Technology Lausanne (EPFL), topo.epfl.ch, Tel:+41(0)21 693 27 55 for providing data, which motivated the research and development of this package. 
#' 
#' @author
#' James Balamuta \email{balamut2@@illinois.edu},
#' Stephane Guerrier \email{stephane@@illinois.edu},
#' Roberto Molinari \email{roberto.molinari@@unige.ch},
#' Wenchao Yang \email{wyang40@@illinois.edu}
#'
#' Stephane Guerrier \email{stephane@@illinois.edu}
#' @name gmwm-package
#' @docType package
#' @useDynLib gmwm
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices gray.colors hcl
#' @importFrom graphics lines plot
#' @importFrom methods is
#' @importFrom stats arima predict ts as.ts
#' @importFrom utils install.packages tail head packageDescription compareVersion
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom grid textGrob gpar
#' @import ggplot2 reshape2
#' @exportPattern ^[[:alpha:]]+
NULL
