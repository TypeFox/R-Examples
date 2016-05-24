#' @name chicks_data
#' @title Data of growth of chicks
#' @description This dataset content dynamic measurements of growth of chicks for different individuals and different strains.
#' The data comes from a selection experiment on chicken initiated by F. Ricard Research Station on Poultry of INRA Nouzilly.
#' The selection focuses on weight at 8 and 36 weeks and allowed to differentiate
#' the following five strains:
#'
#' strain 1 : X-+ (low at 8, but high at 36 weeks)
#'
#' strain 2 : X+- (high at 8, but low at 36 weeks)
#'
#' strain 3 : X++ (high weight at both ages)
#'
#' strain 4 : X- (low weight for both ages)
#'
#' strain 5 : X0 (control).
#'
#'  This is a sub-sample of 50 females in the last generation of selection with weight data (in g) at 12 different ages   to complete measurement (0, 4, 6, 8, 12, 16, 20, 24; 28, 32, 36 and 40 weeks).
#' @docType data
#' @usage carcass_data
#' @format a \code{RangedData} instance, 1 row : strain ; id_animal ; time (day) ; liveweight (g).
#' @source Duval M., Robert-Granie C., Foulley J.-L. (2009) Estimation of heterogeneous variances in non linear mixed models via the SAEM algorithm with applications to growth curves in poultry. Journal de la Societe Francaise de Statistique, 150,65-83
#'
#' Donnet S., Foulley J.-L., Samson A. (2010) Bayesian analysis of growth curves using mixed models defined by stochastic differential equations. Biometrics 66, 733-741
#'
#' Jaffrezic F., Meza C., Foulley .J.-L., Lavielle M. (2006) The SAEM algorithm for the analysis of non linear traits in genetic studies. Genetics, Selection, Evolution,38, 583-
#'
#' This dataset was used in a training session Biobayes (France, 2011) training session.
#'
#' Albert I., Ancelet S., David O., Denis J.B., Makowski D., Parent E., Soubeyrand S. (2012) Methodes statistiques bayesiennes. Bases theoriques et applications en alimentation, environnement et genetique. FormaScience. Ecole-chercheurs INRA.
NULL
