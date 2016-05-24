##' Quality control and reliability
##'
##' This package allows to generate Shewhart-type charts and to obtain 
##' numerical results of interest to the quality control of a process 
##' (involving  continuous, attribute or count data).
##' This package provides basic functionality for univariable and multivariable 
##' quality control analysis, including: xbar, xbar-one, S, R, n, np, c, g, ewna, cusum, 
##' mewna, mcusum  and T2 charts. 
##'
##' @name qcr
##' @aliases qcr
##' @docType package
##' @title Quality Control and Reability
##' @import qcc
NULL

##' @title Vickers hardness data
##' @description A known chemical company is developing a patent for 
##' a new variant of artificial stone composed mostly of quartz ( 93wt %) 
##' and polyester resin . This company is launching a pilot plant where 
##' it begins to produce plates of this material to industry scale. In order 
##' to measure the degree of product homogeneity, 50 samples were taken, 
##' performed 5 measurements per plate corresponding to different areas 
##' of artificial stone Vickers hardness
##'
##' @name plates
##' @docType data
##' @format A data frame with 250 observations on the following 2 variables:
##' \describe{
##'   \item{hardness}{Vickers hardness corresponding to different
##'    areas of artificial stone}
##'   \item{sample}{sample id}
##' }
##' @keywords datasets
##' @examples
##' 
##' data(plates)
##' attach(plates)
##' summary(plates)
##' plot(hardness, type="b")
##' detach(plates)
NULL

##' @title Level of presion data
##' @description A shipyard of recreational boats manufacturing, 
##' intended to optimize and control the mechanical properties hull yacht models. 
##' This has made a study in which the modulus of elasticity tensile strength of the epoxy resin 
##' (polymer) used, after applying different curing pressures measured: 0.1 y 10 MPa.  
##' 60 subsamples composed of three measurements taken on the same day are taken.
##'
##' @name presion
##' @docType data
##' @format A data frame with 180 observations on the following 3 variables:
##' \describe{
##'   \item{presion}{presion level}
##'   \item{sample}{sample id}
##'   \item{measur}{pressures measured: 0.1 y 10 MPa}
##' }
##' @keywords datasets
##' @examples
##' 
##' data(presion)
##' attach(presion)
##' summary(presion)
##' plot(presion$presion, type="b")
##' detach(presion)
NULL

##' @title The performance of the counters data
##' @description A water company from A Corunia wants to control the performance 
##' of the counters installed throughout the city.
##' 60 subsamples are taken each one composed by 3 measurements made 
##' by the counters of the same antiquity (10 years) and caliber, 
##' in a period of 5 years. Taking into account that there are two brands 
##' or providers of counters
##'
##' @name counters
##' @docType data
##' @format A data frame with 180 observations on the following 3 variables:
##' \describe{
##'   \item{error}{the measurement error of the counters (Error:  (Real Volume - Measured Volume)/Real Volume)}
##'   \item{sample}{sample id}
##'   \item{brand}{brands of providers of counters}
##' }
##' @keywords datasets
##' @examples
##' 
##' data(counters)
##' attach(counters)
##' summary(counters)
##' plot(error, type="b")
##' detach(counters)
NULL


##' @title Level of employment data
##' @description A Spanish-Argentina hotel company wants to control 
##' the level of employment in their establishments. 
##' For this it is going to make a continuous control 
##' that measures the amount of occupants in terms of percentage. 
##' 48 sub samples are taken of six hotels belonging to two countries
##'
##' @name employment
##' @docType data
##' @format A data frame with 288 observations on the following 3 variables:
##' \describe{
##'   \item{occupantion}{the amount of occupants in terms of percentage}
##'   \item{sample}{sample id}
##'   \item{hemisphere}{Hemisphere}
##' }
##' @keywords datasets
##' @examples
##' 
##' data(employment)
##' attach(employment)
##' summary(employment)
##' boxplot(occupantion ~ hemisphere)
##' plot(occupantion, type="b")
##' detach(employment)
NULL

##' @title Oxidation Onset Temperature
##' @description This database contains information on the level of purity 
##' of each batch of Picual varities. Then we have the types of oils 
##' by measuring the Oxidation Onset Temperature. 
##' We have 50 subsamples of oils with their temperature to oxide.
##'
##' @name oxidation
##' @docType data
##' @format A data frame with 250 observations on the following 2 variables:
##' \describe{
##'   \item{OOT}{That is a quantitative variable that controls the quality of the oil.}
##'   \item{sample}{sample id}
##' }
##' @keywords datasets
##' @examples
##' 
##' data(oxidation)
##' attach(oxidation)
##' summary(oxidation)
##' plot(OOT, type="b")
##' detach(oxidation)
NULL


##' @title Circuit boards data
##' @description Number of nonconformities observed in 26 successive samples of 100 printed
##' circuit boards. Sample 6 and 20 are outside the control limits. Sample 6
##' was examined by a new inspector and he did not recognize several type of
##' nonconformities that could have been present. Furthermore, the unusually
##' large number of nonconformities in sample 20 resulted from a temperature
##' control problem in the wave soldering machine, which was subsequentely
##' repaired. The last 20 samples are further samples collected on inspection
##' units (each formed by 100 boards).
##' 
##' @name circuit
##' @docType data
##' @format A data frame with 46 observations on the following 4 variables:
##' \describe{
##'   \item{x}{number of defectives in 100 printed circuit boards (inspection unit)}
##'   \item{sample}{sample ID}
##'   \item{size}{sample size}
##'   \item{trial}{trial sample indicator (TRUE/FALSE)}
##' }
##' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
##' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 173--175
##' @keywords datasets
##' @examples
##' 
##' data(circuit)
##' attach(circuit)
##' summary(circuit)
##' boxplot(x ~ trial)
##' plot(x, type="b")
##' detach(circuit)
NULL


##' @title Orange juice data
##' @description Frozen orange juice concentrate is packed in 6-oz cardboard cans. These
##' cans are formed on a machine by spinning them from cardboard stock and
##' attaching a metal bottom panel. A can is then inspected to determine
##' whether, when filled, the liquid could possible leak either on the side
##' seam or around the bottom joint. If this occurs a can is considered
##' nonconforming. The data were collected as 30 samples of 50 cans each at
##' half-hour intervals over a three-shift period in which the machine was in
##' continuous operation. From sample 15 used a new bacth of cardboard stock
##' was punt into production. Sample 23 was obtained when an inexperienced
##' operator was temporarily assigned to the machine. After the first 30
##' samples, a machine adjustment was made. Then further 24 samples were taken
##' from the process.
##'  
##' @name orangejuice
##' @docType data
##' @format A data frame with 54 observations on the following 4 variables:
##' \describe{ 
##' \item{sample}{sample id}
##' \item{D}{number of defectives}
##' \item{size}{sample sizes} 
##' \item{trial}{trial samples (TRUE/FALSE)}
##'  }
##' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
##' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 152--155.
##' @keywords datasets
##' @examples
##' 
##' data(orangejuice)
##' orangejuice$d <- orangejuice$D/orangejuice$size
##' attach(orangejuice)
##' summary(orangejuice)
##' boxplot(d ~ trial)
##' mark <- ifelse(trial, 1, 2)
##' plot(sample, d, type="b", col=mark, pch=mark)
NULL



##' @title Personal computer manufacturer data
##' @description A personal computer manufacturer counts the number of nonconformities per
##' unit on the final assembly line. He collects data on 20 samples of 5
##' computers each.
##' 
##' @name pcmanufact
##' @docType data
##' @format A data frame with 10 observations on the following 2 variables.
##' \describe{
##' \item{x}{number of nonconformities (inspection units)}
##' \item{sample}{sample ID}
##' \item{size}{number of computers inspected}
##' } 
##' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
##' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 181--182
##' @keywords datasets
##' @examples
##' 
##' data(pcmanufact)
##' summary(pcmanufact)
##' plot(pcmanufact$x/pcmanufact$size, type="b")
NULL

##' @title Piston rings data
##' @description Piston rings for an automotive engine are produced by a forging process.
##' The inside diameter of the rings manufactured by the process is measured on
##' 25 samples, each of size 5, drawn from a process being considered 'in
##' control'.
##' 
##' @name pistonrings
##' @docType data
##' @format A data frame with 200 observations on the following 3 variables.
##' \describe{ 
##' \item{diameter}{a numeric vector}
##' \item{sample}{sample ID}
##' \item{trial}{trial sample indicator (TRUE/FALSE)}
##'  }
##' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
##' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 206--213
##' @keywords datasets
##' @examples
##' 
##' data(pistonrings)
##' attach(pistonrings)
##' summary(pistonrings)
##' boxplot(diameter ~ sample)
##' plot(sample, diameter, cex=0.7)
##' lines(tapply(diameter,sample,mean))
##' detach(pistonrings)
NULL

##' @title Target archery dataset in the ranking round (used as Phase I)
##' @description It consists in an stage in which the archer shoots 72 arrows in 12 ends of six arrows. 
##' The information is given in x and y coordinates.
##' 
##' @name archery1
##' @docType data
##' @format An array of (24 x 2 x 3).
##' \describe{ 
##' \item{x-coordinate}{x-coordinate}
##' \item{y-coordinate}{y-coordinate}
##' }
##' @keywords datasets
##' @examples
##'
##' data(archery1)
##' str(archery1) ; plot(archery1)
NULL


##' @title Dowel pin dataset
##' @description Diameter and length of a manufacturing process of a dowel pin##' 
##' 
##' @name dowel1
##' @docType data
##' @format A data frame with 40 observations on the following 2 variables.
##' \describe{ 
##' \item{diameter}{a numeric vector}
##' \item{length}{a numeric vector}
##' }
##' @keywords datasets
##' @examples
##'
##' data(dowel1)
##' str(dowel1) ; plot(dowel1)
NULL
