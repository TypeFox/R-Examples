#############################################################################
#
#  This file is a part of the R package "RoughSets".
#
#  Author: Lala Septem Riza and Andrzej Janusz
#  Supervisors: Chris Cornelis, Francisco Herrera, Dominik Slezak and Jose Manuel Benitez
#  Copyright (c):
#       DiCITS Lab, Sci2s group, DECSAI, University of Granada and
#       Institute of Mathematics, University of Warsaw
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' Several datasets have been embedded in this package to be used as decision table of examples.
#' They can be accessed by typing \code{data(RoughSetData)}. The following is a description of each 
#' datasets.
#'
#' \bold{The hiring dataset}
#' 
#' It is simple data taken from (Komorowski et al, 1999) where all the attributes have nominal values. 
#' It consists of eight objects with four conditional attributes and one decision attribute. 
#' The detailed description of each attribute is as follows:
#' \itemize{
#' \item Diploma: it has the following values: {"MBA", "MSc", "MCE"}.
#' \item Exprience: it has the following values: {"High", "Low", "Medium"}.
#' \item French: it has the following values: {"Yes", "No"}.
#' \item Reference: it has the following values: {"Excellent", "Good", "Neutral"}.
#' \item Decision: it is a decision attribute that contains the following values: {"Accept", "Reject"}.
#' }
#'
#' \bold{The housing dataset}
#'
#' This data was taken from the Boston housing dataset located at the UCI Machine Learning repository, 
#' available at http://www.ics.uci.edu. It was first created by 
#' (Harrison and Rubinfeld, 1978). It contains 506 objects with 13 conditional attributes and one decision attribute.
#' Furthermore, it should be noted that the housing dataset is a regression dataset which means that 
#' the decision attribute has continuous values. The conditional attributes contain both continuous and nominal attributes.
#' The following is a description of each attribute:
#' \itemize{
#' \item CRIM: it is a continuous attribute that expresses per capita crime rate by town. 
#'             It has values in: [0.0062, 88.9762].
#' \item ZN: it is a continuous attribute that represents the proportion of residential land 
#'             zoned for lots over 25,000 sq.ft. It has values in: [0, 100].
#' \item INDUS: it is a continuous attribute that shows the proportion of non-retail business acres per town.
#'              It has values in: [0.46, 27.74].
#' \item CHAS: it is a nominal attribute that represents Charles River dummy variable. 
#'              It has two values which are 1 if tract bounds river and 0 otherwise.
#' \item NOX: it is a continuous attribute that shows the nitric oxides concentration (parts per 10 million).
#'            It has values in: [0.385, 0.871].
#' \item RM: it is a continuous attribute that explains the average number of rooms per dwelling.
#'           It has values in: [3.561, 8.78].
#' \item AGE: it is a continuous attribute that expresses proportion of owner-occupied 
#'           units built prior to 1940. It has values in: [2.9, 100].
#' \item DIS: it is a continuous attribute that shows weighted distances to five Boston employment centres.
#'            It has values in: [1.1296, 12.1265].
#' \item RAD: it is a nominal attribute that shows the index of accessibility to radial highways.
#'            it has the integer value from 1 to 24.
#' \item TAX: it is a continuous attribute that shows the full-value property-tax rate per $10,000.
#'            It has values in: [187, 711].
#' \item PTRATIO: it is a continuous attribute that shows the pupil-teacher ratio by town.
#'                It has values in: [12.6, 22].
#' \item B: it is a continuous attribute that can be expressed by 1000(Bk - 0.63)^2 
#'          where Bk is the proportion of blacks by town. It has values in: [0.32, 396.9].
#' \item LSTAT: it is a continuous attribute that illustrates the percentage of lower status of the population.
#'          It has values in: [1.73, 37.97].
#' \item MEDV: it is a continuous attribute that shows the median value of owner-occupied homes in $1000's.
#'             It has values in: [5, 50].
#' }
#'
#' \bold{The wine dataset}
#' 
#' This dataset is a classification dataset introduced first by (Forina, et al) 
#' which is commonly used as benchmark for simulation in the machine learning area. 
#' Additionally, it is available at the KEEL dataset repository (Alcala-Fdez, 2009), available at http://www.keel.es/. 
#' It consists of 178 instances with 13 conditional attributes and one decision attribute 
#' where all conditional attributes have continuous values. The description of 
#' each attribute is as follows:
#' \itemize{
#' \item alcohol: it has a range in: [11, 14.9].
#' \item malid_acid: it has a range in: [0.7, 5.8].
#' \item ash: it has a range in: [1.3, 3.3].
#' \item alcalinity_of_ash: it has a range in: [10.6, 30.0].
#' \item magnesium: it has a range in: [70, 162].
#' \item total_phenols: it has a range in: [0.9, 3.9].
#' \item flavanoids: it has a range in: [0.3 5.1].
#' \item nonflavanoid_phenols: it has a range in: [0.4 3.6].
#' \item proanthocyanins: it has a range in: [0.4 3.6].
#' \item color_intensity: it has a range in: [1.2 13.0].
#' \item hue: it has a range in: [0.4 1.8].
#' \item od: it has a range in: [1.2 4.0].
#' \item proline: it has a range in: [278 1680].
#' \item class: it is nominal decision attribute that has values: {1, 2, 3}.
#' }
#'
#' \bold{The pima dataset}
#' 
#' It was taken from the pima Indians diabetes dataset which is available at the KEEL dataset repository 
#' (Alcala-Fdez, 2009), available at http://www.keel.es/. It was first created by 
#' National Institute of Diabetes and Digestive and Kidney Diseases. It contains 768 objects with 8 continuous conditional
#' attributes. The description of each attribute is as follows:
#' \itemize{
#' \item preg: it represents number of times pregnant and has values in: [1, 17].
#' \item plas: it represents plasma glucose concentration 
#'              a 2 hours in an oral glucose tolerance test and has values in: [0.0, 199.0].
#' \item pres: it represents diastolic blood pressure (mm Hg) and has values in: [0.0, 122.0].
#' \item skin: it represents triceps skin fold thickness (mm) and has values in: [0.0, 99.0].
#' \item insu: it represents 2-hour serum insulin (mu U/ml) and has values in: [0.0, 846.0].
#' \item mass: it represents body mass index (weight in kg/(height in m)^2) and has values in: [0.0, 67.1].
#' \item pedi: it represents diabetes pedigree function and has values in: [0.078, 2.42].
#' \item age: it represents age (years) and has values in: [21.0, 81.0].
#' \item class: it is a decision attribute and has values in: [1, 2].
#' }
#' 
#' @title Data set of the package
#' @name RoughSetData
#' @docType data
#' @references 
#' 
#' M. Forina, E. Leardi, C. Armanino, and S. Lanteri, "PARVUS -
#' An Extendible Package for Data Exploration, Classification and Correlation",
#' Journal of Chemonetrics, vol. 4, no. 2, p. 191 - 193 (1988). 
#'
#' D. Harrison, and D. L. Rubinfeld, "Hedonic Prices and the 
#' Demand for Clean Air", J. Environ. Economics & Management,
#' vol.5, 81-102 (1978).
#' 
#' J. Alcala-Fdez, L. Sanchez, S. Garcia, M. J. del Jesus, S. Ventura, 
#' J. M. Garrell, J. Otero, C. Romero, J. Bacardit, V. M. Rivas, 
#' J. C. Fernandez, and F. Herrera,  
#' "KEEL: A Software Tool to Assess Evolutionary Algorithms to Data Mining Problems", 
#' Soft Computing vol. 13, no. 3, p. 307 - 318 (2009).
#'
#' J. Komorowski, Z. Pawlak, L. Polwski, and A. Skowron,
#' "Rough Sets: A Tutorial", In S. K. Pal and A. Skowron, editors,
#' Rough Fuzzy Hybridization, A New Trend in Decision Making,
#' pp. 3 - 98, Singopore, Springer (1999).
#' 
#' @keywords data
NULL