# Copyright (C) 2012-2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or
#     modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see
#     <http://www.gnu.org/licenses/>.

#' @import MASS
#' @import lhs
#' @import plyr
#' @import BBmisc
#' @import ParamHelpers
#' @import ggplot2
#' @import metafor
#' @import BatchJobs
#' @import BatchExperiments
NULL

if (getRversion() >= '2.15.1') globalVariables(c('tpos', 'tneg', 'cpos',
                                                 'cneg', 'dat.bcg',
                                                 'type', 'bias',
                                                 'confidence',
                                                 'coverage', 'method',
                                                 'mse', 'variance',
                                                 'width', 'lower',
                                                 'upper',
                                                 'heterogeneity',
                                                 '..density..',
                                                 '..count..',
                                                 'intercept', 'via',
                                                 'slope', 'reference',
                                                 'logrisk', 'rlower',
                                                 'rupper', 'size',
                                                 'balance', 'h',
                                                 'muted', 'colour',
                                                 'd_mean',
                                                 'cnfColourPalette',
                                                 'estColourPalette',
                                                 'd_sd'))
