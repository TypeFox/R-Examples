#' Additive and Multiplicative Effects Models for Networks and Relational Data
#' 
#' Analysis of network and relational data using additive and multiplicative
#' effects (AME) models. The basic model includes regression terms,
#' the covariance structure of the social relations model
#' (Warner, Kenny and Stoto (1979), Wong (1982)), and multiplicative
#' factor effects (Hoff(2009)). Four different link functions accommodate
#' different relational data structures, including binary/network data (bin),
#' normal relational data (nrm), ordinal relational data (ord) and data from
#' fixed-rank nomination schemes (frn).  Several of these link functions are
#' discussed in Hoff, Fosdick, Volfovsky and Stovel (2013). Development of this
#' software was supported in part by NICHD grant R01HD067509. 
#' 
#' \tabular{ll}{ Package: \tab amen\cr Type: \tab Package\cr Version: \tab
#' 1.1 \cr Date: \tab 2015-05-31 \cr License: \tab GPL-3 \cr }
#' 
#' @name amen-package
#' @aliases amen-package amen
#' @docType package
#' @author Peter Hoff, Bailey Fosdick, Alex Volfovsky, Yanjun He
#' 
#' Maintainer: Peter Hoff <pdhoff@@uw.edu>
#' @keywords package
#' @examples
#' 
#' 
#' data(YX_frn)
#' fit<-ame(YX_frn$Y,YX_frn$X,burn=5,nscan=5,odens=1,model="frn")
#' 
#' summary(fit)
#' 
#' plot(fit) 
#' 
#' 
NULL




#' binary relational data and covariates
#' 
#' a synthetic dataset that includes binary relational data as well as
#' information on eight covariates
#' 
#' 
#' @name YX_bin
#' @docType data
#' @usage data(YX_bin)
#' @format The format is: List of 2 $ Y: num [1:100, 1:100] NA 0 0 0 0 0 0 0 0
#' 1 ...  $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...  ..- attr(*,
#' "dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr [1:8]
#' "intercept" "rgpa" "rsmoke" "cgpa" ...
#' @keywords datasets
#' @examples
#' 
#' data(YX_bin)
#' gofstats(YX_bin$Y) 
#' 
NULL





#' Censored binary nomination data and covariates
#' 
#' a synthetic dataset that includes relational data where the number of
#' nominations per row is censored at 10, along with information on eight
#' covariates
#' 
#' 
#' @name YX_cbin
#' @docType data
#' @usage data(YX_cbin)
#' @format The format is: List of 2 $ Y: num [1:100, 1:100] NA 0 0 0 1 0 0 0 0
#' 3 ...  $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...  ..- attr(*,
#' "dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr [1:8]
#' "intercept" "rgpa" "rsmoke" "cgpa" ...
#' @keywords datasets
#' @examples
#' 
#' data(YX_cbin)
#' gofstats(YX_cbin$Y) 
#' 
NULL





#' Fixed rank nomination data and covariates
#' 
#' a synthetic dataset that includes fixed rank nomination data as well as
#' information on eight covariates
#' 
#' 
#' @name YX_frn
#' @docType data
#' @usage data(YX_frn)
#' @format The format is: List of 2 $ Y: num [1:100, 1:100] NA 0 0 0 1 0 0 0 0
#' 3 ...  $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...  ..- attr(*,
#' "dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr [1:8]
#' "intercept" "rgpa" "rsmoke" "cgpa" ...
#' @keywords datasets
#' @examples
#' 
#' data(YX_frn)
#' gofstats(YX_frn$Y) 
#' 
NULL





#' normal relational data and covariates
#' 
#' a synthetic dataset that includes continuous (normal) relational data as
#' well as information on eight covariates
#' 
#' 
#' @name YX_nrm
#' @docType data
#' @usage data(YX_nrm)
#' @format The format is: List of 2 $ Y: num [1:100, 1:100] NA -4.05 -0.181
#' -3.053 -1.579 ...  $ X: num [1:100, 1:100, 1:8] 1 1 1 1 1 1 1 1 1 1 ...  ..-
#' attr(*, "dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr [1:8]
#' "intercept" "rgpa" "rsmoke" "cgpa" ...
#' @keywords datasets
#' @examples
#' 
#' data(YX_nrm)
#' gofstats(YX_nrm$Y)
#' 
#' 
NULL





#' ordinal relational data and covariates
#' 
#' a synthetic dataset that includes ordinal relational data as well as
#' information on seven covariates
#' 
#' 
#' @name YX_ord
#' @docType data
#' @usage data(YX_ord)
#' @format The format is: List of 2 $ Y: num [1:100, 1:100] NA 0 3 0 3 1 0 1 1
#' 0 ...  $ X: num [1:100, 1:100, 1:7] 1 1 1 1 1 1 1 1 1 1 ...  ..- attr(*,
#' "dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr [1:7] "rgpa"
#' "rsmoke" "cgpa" "csmoke" ...
#' @keywords datasets
#' @examples
#' 
#' data(YX_ord)
#' gofstats(YX_ord$Y)
#' 
NULL





#' row-specific ordinal relational data and covariates
#' 
#' a synthetic dataset that includes row-specific ordinal relational data as
#' well as information on five covariates
#' 
#' 
#' @name YX_rrl
#' @docType data
#' @usage data(YX_rrl)
#' @format The format is: List of 2 $ Y: num [1:100, 1:100] NA 0 3 0 3 1 0 1 1
#' 0 ...  $ X: num [1:100, 1:100, 1:5] 1 1 1 1 1 1 1 1 1 1 ...  ..- attr(*,
#' "dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr [1:5] "cgpa"
#' "csmoke" "igrade" "ismoke" ...
#' @keywords datasets
#' @examples
#' 
#' data(YX_rrl)
#' gofstats(YX_rrl$Y)
#' 
NULL


#' binary relational data and covariates
#' 
#' a synthetic dataset that includes longitudinal binary relational data
#' as well as information on covariates
#' 
#' 
#' @name YX_bin_long
#' @docType data
#' @usage data(YX_bin_long)
#' @format a list
#' @keywords datasets
#' @examples
#' 
#' data(YX_bin_long)
#' gofstats(YX_bin_long$Y[,,1]) 
#' 
NULL


#' @title Sampson's monastery data
#'
#' @description 
#' Several dyadic variables measured on 18 members of a monastery. 
#'
#' @format 
#' A socioarray whose dimensions represent nominators, nominatees and relations. 
#' Each monk was asked to rank up to three other monks on a variety of positive 
#' and negative relations. A rank of three indicates the "highest" ranking for 
#' a particular relational variable. The relations \code{like_m2} and \code{like_m1}
#' are evaluations of likeing at one and two timepoints previous to when the 
#' other relations were measured. 
#' 
#' @source 
#' \url{http://moreno.ss.uci.edu/data.html#sampson}
#' 
#' @name sampsonmonks
NULL


#' @title Cold War data
#'
#' @description 
#' Positive and negative relations between countries during the cold war
#'
#' @format 
#' A list including the following dyadic and nodal variables:
#' \itemize{
#' \item \code{cc}: a socioarray of ordinal levels of military 
#' cooperation (positive) and conflict (negative), every 5 years; 
#' \item \code{distance}: between-country distance (in thousands of kilometers); 
#' \item \code{gdp}: country gdp in dollars every 5 years; 
#' \item \code{polity}: country polity every 5 years.  
#' }
#' @source 
#' Xun Cao : \url{http://polisci.la.psu.edu/people/xuc11}
#' 
#' @name coldwar
NULL


#' @title Comtrade data
#'
#' @description
#' Eleven years of import and export data between 229 countries. 
#' The data use the SITC Rev. 1 commodity classification, aggregated at the 
#' first level (AG1). 
#'  
#' @format A list consisting of a socioarray \code{Trade} and a vector 
#' \code{dollars2010} of inflation rates. The socioarray gives 
#' yearly trade volume (exports and imports) 
#' in dollars for 10 different commodity classes 
#' for eleven years between 229 countries. This gives a five-way 
#' array. The first index is the reporting country, so 
#' \code{Trade[i,j,t,k,1]} is what \code{i} reports for exports to 
#' \code{j}, but in general this is not the same as 
#' \code{Trade[j,i,t,k,2]}, what \code{j} reports as importing from \code{i}. 
#' 
#' @source \url{http://comtrade.un.org/}, \url{http://www.measuringworth.com/}
#' 
#' @name comtrade
NULL


#' @title Lazega's law firm data
#'
#' @description 
#' Several nodal and dyadic variables measured on 71 attorneys in a law firm. 
#'
#' @format 
#' A list consisting of a socioarray \code{Y} and a nodal attribute matrix \code{X}. 
#' 
#' The dyadic variables in \code{Y} include three binary networks: advice, friendship
#' and co-worker status. 
#' 
#' The categorical nodal attributes in \code{X} are coded as follows: 
#' \itemize{
#' \item status (1=partner, 2=associate)
#' \item office (1=Boston, 2=Hartford, 3=Providence) 
#' \item practice (1=litigation, 2=corporate) 
#' \item law school (1=Harvard or Yale, 2=UConn, 3=other)  
#'  }
#' \code{seniority} and \code{age} are given in years, and \code{female} is 
#' a binary indicator. 
#' 
#' @source 
#' \url{http://moreno.ss.uci.edu/data.html#lazega}
#' 
#' @name lazegalaw
NULL


#' AddHealth community 3 data
#'
#' A valued sociomatrix (Y) and matrix of nodal attributes (X) for 
#' students in community 3 of the AddHealth study. 
#' \itemize{
#' \item Y: A sociomatrix in which the value of the edge corresponds to an ad-hoc measure of intensity of the relation. Note that students were only allowed to nominate up to 5 male friends and 5 female friends. 
#' \item X: Matrix of students attributes, including sex, race (1=white, 2=black, 3=hispanic, 4=asian, 5=mixed/other) and grade. 
#' } 
#' See \url{http://moreno.ss.uci.edu/data.html#adhealth} for more details. 
#' @docType data
#' @keywords datasets
#' @format list 
#' @name addhealthc3 
#' @usage data(addhealthc3)
NULL

#' AddHealth community 9 data
#'
#' A valued sociomatrix (Y) and matrix of nodal attributes (X) for 
#' students in community 9 of the AddHealth study. 
#' \itemize{
#' \item Y: A sociomatrix in which the value of the edge corresponds to an ad-hoc measure of intensity of the relation. Note that students were only allowed to nominate up to 5 male friends and 5 female friends. 
#' \item X: Matrix of students attributes, including sex, race (1=white, 2=black, 3=hispanic, 4=asian, 5=mixed/other) and grade. 
#' } 
#' See \url{http://moreno.ss.uci.edu/data.html#adhealth} for more details.
#' @docType data
#' @keywords datasets
#' @format list 
#' @name addhealthc9 
#' @usage data(addhealthc9)
NULL



#' @title International relations in the 90s
#'
#' @description 
#' A relational dataset recording 
#' a variety of nodal and dyadic variables on countries in the 1990s,
#' including information on conflicts, trade and other variables. 
#' Except for the conflict variable, the variables are averages
#' across the decade. 
#'
#' @format 
#' A list consisting of a socioarray \code{dyadvars} of 
#' dyadic variables and matrix  \code{nodevars} of nodal variables. 
#' The dyadic variables include
#' \itemize{
#' \item total number of conflicts;
#' \item exports (in billions of dollars);
#' \item distance (in thousands of kilometers);
#' \item number of shared IGOs (averages across the years);
#' \item polity interaction. 
#' } 
#' The nodal variables include
#' \itemize{
#' \item population (in millions); 
#' \item gdp (in billions of dollars); 
#' \item polity 
#' }
#' 
#' @source Michael Ward.
#' 
#' @name IR90s
NULL



#' @title Dutch college data
#'
#' @description
#' Longitudinal relational measurements and nodal characteristics 
#' of Dutch college students, described in 
#' van de Bunt, van Duijn, and Snijders (1999). 
#' The time interval between the first four measurements was 
#' three weeks, whereas the interval between the last three 
#' was six weeks. 
#' 
#' @format A list consisting of a socioarray \code{Y} and a matrix 
#' \code{X}  of static nodal attributes. The relational 
#' measurements range from -1 to 4, indicating the following:
#' \itemize{
#' \item -1 a troubled or negative relationship
#' \item  0 don't know
#' \item  1 neutral relationship
#' \item  2 friendly
#' \item  3 friendship
#' \item  4 best friends
#' }
#' 
#' @source \url{http://moreno.ss.uci.edu/data.html#vdb}
#' 
#' @name dutchcollege
NULL




#' @title Sheep dominance data
#'
#' @description 
#' Number of dominance encounters between 28 female bighorn sheep. 
#' Cell (i,j) records the number of times sheep i dominated sheep j. 
#' From Hass (1991). 
#'  
#' @format 
#' A list consisting of the following:
#' \itemize{
#' \item \code{dom}: a directed socioarray recording the number of 
#' dominance encounters. 
#' \item \code{age}: the age of each sheep in years. 
#' }
#' 
#' @source 
#' \url{http://moreno.ss.uci.edu/data.html#sheep}
#' @name sheep
NULL


