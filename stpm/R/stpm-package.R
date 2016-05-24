#'Stochastic Process Model for Analysis of Longitudinal and Time-to-Event Outcomes
#'
#'Utilities to estimate parameters of the models with survival functions induced by 
#'stochastic covariates. Miscellaneous functions for data preparation and simulation 
#'are also provided. For more information, see: "Stochastic model for analysis of 
#'longitudinal data on aging and mortality" by Yashin A. et al, 2007, Mathematical 
#'Biosciences, 208(2), 538-551 <DOI:10.1016/j.mbs.2006.11.006>.
#'@references Yashin, A. et al (2007), Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.
#'@references Akushevich I., Kulminski A. and Manton K. (2005). Life tables with covariates: Dynamic model 
#'for Nonlinear Analysis of Longitudinal Data. Mathematical Popu-lation Studies, 12(2), pp.: 51-80.
#'<DOI: 10.1080/08898480590932296>.
#'@references Yashin, A. et al (2007), Health decline, aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#'@author I. Y. Zhbannikov, I. V. Akushevich, K. G. Arbeev, A. I. Yashin.
#'@maintainer Ilya Y. Zhbannikov <ilya.zhbannikov@duke.edu>
#'@doctype package
#'@name stpm
#'@aliases stpm stpm-package
#'@keywords Stochastic process model, allostatic load, quadratic hazard
#'@version 1.2.1
#'@date 2016-04-10
#'longitudinal studies, time-to-event
#'@examples \dontrun{ 
#'library(stpm)
#'#Prepare data for optimization
#'data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), 
#'					   y=system.file("data","vitstat.csv",package="stpm"))
#'#Parameters estimation (default model: discrete-time):
#'p.discr.model <- spm(data)
#'p.discr.model
#'# Continuous-time model:
#'p.cont.model <- spm(data, model="continuous")
#'p.cont.model
#'#Model with time-dependent coefficients:
#'data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), 
#'					   y=system.file("data","vitstat.csv",package="stpm"), 
#'					   covariates="BMI")
#'p.td.model <- spm(data, model="time-dependent")
#'p.td.model
#'}
NULL