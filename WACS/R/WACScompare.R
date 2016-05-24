###################################################################
#
# This function is part of WACSgen 1.0
# Copyright © 2013,2014,2015, D. Allard, BioSP,
# and Ronan Trepos MIA-T, INRA
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. http://www.gnu.org
#
###################################################################

#' Performs comparisons between two WACS data structures, or between two WACS simulation series
#' 
#' The comparison is based on different types of 
#' statistics computed on WACSdata1 and WACSdata2, or WACSsim1 and WACSsim2
#' 
#' @export
#' 
#' @param what      Type of validation. Possible choices are:
#'  \tabular{ll}{
#'  \env{="Sim"} \tab Compares a simulation run to data \cr
#'  \env{="Rain"}\tab qq-plots of rainfall, per season \cr
#'  \env{="MeanSd"}\tab Compares monthly mean and standard deviations \cr
#'  \env{="BiVar"}\tab Compares bivariate correlations \cr
#'  \env{="CorTemp"}\tab Compares temporal correlations \cr
#'  \env{="SumBase"}\tab Compares sums above a threshold \cr
#'  \env{="Persistence"}\tab Compares persistence of a variable above (or below) a threshold \cr
#'  }
#' @param wacs1  Either WACS data obtained when calling \link{WACSdata}, or WACS simulations obtained
#'  when calling \link{WACSsimul}.
#'  
#' @param wacs2 Either WACS data obtained when calling \link{WACSdata}, or WACS simulations obtained
#'  when calling \link{WACSsimul}. Must be of the same class as \code{wacs1}
#'  
#' @param wacspar   WACS parameters estimated when calling \link{WACSestim} on \code{wacs1}
#'  
#' @param varname   Variable on which the validation is performed
#'  
#' @param varname2  Second variable on which validation is performed (only needed if \code{what=BiVar})
#'  
#' @param base      Threshold used for "SumBase" and "Persistence"
#'  
#' @param above     Boolean value used for "Persistence": 
#'    \code{TRUE} if data is considered above threshold; 
#'    \code{FALSE} otherwise
#'  
#' @param months    Months to which the analysis should be restricted (only for "SumBase" and "Persistence")
#'
#' @return A list containing all information needed for plots; contains also the type of validation, as a class
#'  
#' @examples 
#' \dontrun{
#'   ## Simple example
#'   data(ClimateSeries)
#'   ThisData = WACSdata(ClimateSeries,from="1995-01-01",to=1999-12-31")
#'   ThisPar  = WACSestim(ThisData)
#'   ThatData = WACSdata(ClimateSeries,from="2002-01-01",to="1999-12-31")
#'   Comp = WACScompare(what="Sim", wacs1=ThisData, wacspar=ThisPar,
#'                      wacs2=ThatData, varname="tmin")
#'   WACSplot(Comp)
#'   Comp = WACScompare(what="MeanSd",wacs1=ThisData, wacspar=ThisPar,
#'                      wacs2=ThatData, varname="RG")
#'   WACSplot(Comp)
#'   Comp = WACScompare(what="SumBase", wacs1=ThisData, wacspar=ThisPar,
#'                      wacs2=ThatData, varname="tmoy", base=5, months=2:5)
#'   WACSplot(Comp)
#'   Comp = WACScompare(what="Persistence",wacs1=ThisData, wacspar=ThisPar,
#'                      wacs2=ThatData, varname="tmin", base=0, above=FALSE)
#'   WACSplot(Comp)
#' }

#'  
#' @note
#' 
#' \code{wacs1} and \code{wacs2} must be of the same \code{WACS} class. 
#' We must either have 
#' 
#' \code{class(wacs1)=class(wacs2)=class(WACSdata)}, 
#' 
#' or \code{class(wacs1)=class(wacs2)=class(WACSsimul)}.
#' 
#' If \code{what="sim"}, data and simulations are displayed as a function of the day of the year, from 1 to 365. 
#' Smoothed versions daily average and daily envelopes (defined by average +/- 2. standard deviations) are also displayed. 
#' 
#' 
#' If \code{what="rain"}, qq-plots and superimposition of histograms and models of rain are produced for each season.
#' 
#' 
#' If \code{what="MeanSd"}, boxplots of monthly means and monthly standard deviations are compared. 
#' The median value of the monthly mean, resp. monthly standard deviation, of the data are displayed 
#' on top of the boxplots computed on the simulations.
#' 
#' 
#' If \code{what="BiVar"}, boxplots of monthly correlations coefficients between \code{varname} and \code{varname2}
#' are compared. The median value of the correlation coefficient computed on the data is displayed 
#' on top of the boxplots computed on the simulations.
#' 
#' 
#' If \code{what="CorTemp"}, boxplots of monthly 1-day auto correlation are compared. The median value of the 
#' auto-correlation coefficient computed on the data is displayed 
#' on top of the boxplots computed on the simulations.
#' 
#' 
#' If \code{what="SumBase"}, boxplots of the sum of the variable \code{varname} above a given threshold, \code{base},
#' is computed during the months provided in the variable \code{months}.
#' 
#' 
#' If \code{what="Persistence"}, histograms of consecutive days of the variable \code{varname} above (or below) a given
#' threshold,  \code{base}, are compared. If \code{above=TRUE}, consecutive days above the threshold are computed, 
#' whereas days below the threshold are computed if \code{above=FALSE}. Months can be selected 
#' with the variable \code{months}.
#' 
#' 

WACScompare = function(what=what,
                     wacs1 = wacs1, 
                     wacspar   = wacspar,
                     wacs2 = wacs2, 
                     varname=varname, 
                     varname2=NULL, 
                     base = 0, 
                     above=T, 
                     months=1:12)
{
  # Some checking
  
  if (!(what %in% c("Sim","Rain","MeanSd","BiVar","CorTemp","SumBase","Persistence"))){
    stop ("[WACScompare: 'what' should belong to one of the following: 
          'Sim','Rain','MeanSd','BiVar','CorTemp','SumBase','Persistence']")
  } 
  
  if(class(wacs1)!=class(wacs2)){
    stop ("[WACScompare] Data should be of same class. Either 'WACSdata', as generated by calling WACSdata or 'WACSsimu' 
           as generated by calling WACSsimu")
  }
  
  if( (class(wacs1) != "WACSdata") & (class(wacs1) != "WACSsimul") ){
    stop("[WACScompare: Input should be of class 'WACSdata', as generated by calling WACSdata, or of class 'WACSsimu' 
           as generated by calling WACSsimu")
  }
  
  
  
  if (what != "Rain"){
    if (!(varname %in% wacspar$mapping$wacs_names) & !(varname=="tmoy")){
      stop ("[WACScompare: varname must be one of the variable name of WACSdata/WACSsimu]")
    }
    
    # Creating 'tmax' if Trange = TRUE
    if( ((varname=="tmax") && wacspar$Trange)  ||  ((varname=="tmoy") && wacspar$Trange)  ){
      tmax = wacs1$data$tmin + wacs1$data$trange
      wacs1$data = cbind(wacs1$data,tmax)
      tmax = wacs2$data$tmin + wacs2$data$trange
      wacs2$data = cbind(wacs2$data,tmax)
    }
  }
  wacsdata = NULL;
  wacssimul = NULL;
  
  if (class(wacs1) == "WACSdata"){
    wacsdata = wacs1
    wacssimul = list(sim=wacs2$data)
    class(wacssimul) = "WACSsimul"
  }
  if (class(wacs1) == "WACSsimul"){
    wacsdata = list(data = wacs1$sim,mapping = wacspar$mapping, bounds=wacspar$bounds, 
                     seasons= wacspar$seasons, Trange = wacspar$Trange)
    wacsdata$data = wacsdata$data[-(which(wacsdata$data$month == 2 & wacsdata$data$day == 29)), ]
    class(wacsdata)="WACSdata"
    wacssimul=wacs2
  }

  res = WACSvalid(what=what,wacsdata=wacsdata,wacspar=wacspar,wacssimul=wacssimul,
                  varname=varname,varname2=varname2,base=base,above=above,months=months)
  if (class(wacs1) == "WACSdata"){
    res$labels = c("Observed 1","Observed 2")
  }
  if (class(wacs1) == "WACSsimul"){
    res$labels = c("Simulations 1","Simulations 2")
  }

  return(res)
}
