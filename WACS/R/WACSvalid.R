###################################################################
#
# This function is part of WACSgen 1.0
# Copyright © 2013,2014,2015, D. Allard, BioSP,
# and Ronan Trépos MIA-T, INRA
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

#' Performs validations of WACS simulations 
#' 
#' The validation is based on different types of 
#' statistics computed on  WACS data, WACS parameters and WACS simulations.
#'
#' @export
#'
#' @param what      Type of validation. Possible choices are:
#'  \tabular{ll}{
#'  \env{="Sim"} \tab Compares a simulation run to data \cr
#'  \env{="Rain"}\tab qq-plots of rainfall, per season \cr
#'  \env{="MeanSd"}\tab Compares monthly mean and standard deviations \cr
#'  \env{="BiVar"}\tab Compares monthly variate correlations \cr
#'  \env{="CorTemp"}\tab Compares monthly temporal correlations \cr
#'  \env{="SumBase"}\tab Compares sums above a threshold \cr
#'  \env{="Persistence"}\tab Compares persistence of a variable above (or below) a threshold \cr
#'  }
#' @param wacsdata  WACS data obtained when calling \link{WACSdata}
#'  
#' @param wacspar   WACS parameters estimated when calling \link{WACSestim}
#'  
#' @param wacssimul WACS simulation obtained when calling \link{WACSsimul}
#'  
#' @param varname   Variable on which the validation is performed
#'  
#' @param varname2  Second variable on which validation is performed (only needed if \code{what=BiVar})
#' 
#' @param base      Threshold used for "SumBase" and "Persistence"
#'  
#' @param above     Boolean value used for "Persistence": 
#'    TRUE if data is considered above threshold; 
#'    FALSE otherwise
#'    
#' @param months    Months to which the analysis should be restricted (only for "SumBase" and "Persistence")
#'
#' @return A list containing all information needed for plots; contains also the type of validation, as a class
#'  
#' @examples 
#' \dontrun{
#'   ## Simple example
#'   data(ClimateSeries)
#'   ThisData = WACSdata(ClimateSeries,from="1995-01-01",to="1999-12-31")
#'   ThisPar  = WACSestim(ThisData)
#'   ThisSim  = WACSsimul(ThisPar,from="1995-01-01",to="1999-12-31")
#'   Val1 = WACSvalid(what="Sim",wacsdata = ThisData, 
#'                    wacspar = ThisPar, wacssimul = ThisSim, varname="tmin")
#'   Val2 = WACSvalid(what="MeanSd",wacsdata = ThisData, 
#'                    wacssimul = ThisSim, varname="RG")
#'   Val3 = WACSvalid(what="SumBase", wacsdata = ThisData,  
#'                    wacssimul = ThisSim, varname="tmoy", base=5, month=2:5)
#'   Val4 = WACSvalid(what="Persistence",wacsdata = ThisData,  
#'                    wacssimul = ThisSim, varname="tmin", base=0, above=FALSE)
#' }
#' @note
#' 
#' If \code{what=sim}, data and simulations are displayed as a function of the day of the year, from 1 to 365. 
#' Smoothed versions of daily average and daily envelopes (defined by average +/- 2. standard deviations) are also displayed. 
#' 
#' 
#' If \code{what=rain}, qq-plots and superimposition of histograms and models of rain are produced for each season.
#' 
#' 
#' If \code{what=MeanSd}, boxplots of monthly means and monthly standard deviations are compared. 
#' The median value of the monthly mean, resp. monthly standard deviation, of the data are displayed 
#' on top of the boxplots computed on the simulations.
#' 
#' 
#' #' If \code{what=BiVar}, boxplots of monthly correlations coefficients between \code{varname} and \code{varname2}
#' are compared. The median value of the correlation coefficient computed on the data is displayed 
#' on top of the boxplots computed on the simulations.
#' 
#' 
#' If \code{what=CorTemp}, boxplots of monthly 1-day auto correlation are compared. The median value of the 
#' auto-correlation coefficient computed on the data is displayed 
#' on top of the boxplots computed on the simulations.
#' 
#' 
#' If \code{what=SumBase}, boxplots of the sum of the variable \code{varname} above a given threshold, \code{base},
#' is computed during the months provided in the variable \code{months}.
#' 
#' 
#' If \code{what=Persistence}, histograms of consecutive days of the variable \code{varname} above (or below) a given
#' threshold,  \code{base}, are compared. If \code{above=TRUE}, consecutive days above the threshold are computed, 
#' whereas days below the threshold are computed if \code{above=FALSE}. Months can be selected 
#' with the variable \code{months}.
#' 
#' 



WACSvalid = function(what="Sim",
                     wacsdata=NULL, 
                     wacspar=NULL, 
                     wacssimul=NULL, 
                     varname=NULL, 
                     varname2=NULL, 
                     base = 0, 
                     above=TRUE, 
                     months=1:12)
{
  # Some checking
  
  if (!(what %in% c("Sim","Rain","MeanSd","BiVar","CorTemp","SumBase","Persistence"))){
    stop ("[WACSvalid: 'what' should belong to one of the following: 
          'Sim','Rain','MeanSd','BiVar','CorTemp','SumBase','Persistence']")
  } 
  if(class(wacsdata)!="WACSdata"){
    stop ("[WACSvalid] Data should be of class 'WACSdata', as generated by calling WACSdata")
  }
  if(!is.null(wacssimul) && class(wacssimul)!="WACSsimul"){
    stop ("[WACSvalid] Simulation should be of class 'WACSsimul', as generated by calling WACSsimul")
  } 
  if (what != "Rain"){
    if (!(varname %in% wacsdata$mapping$wacs_names) & !(varname=="tmoy")){
      stop ("[WACSvalid: varname must be one of the variable name of WACSdata]")
    }
    
    # Creating 'tmax' if Trange = TRUE
    if( ((varname=="tmax") && wacsdata$Trange)  ||  ((varname=="tmoy") && wacsdata$Trange)  ){
      tmax = wacsdata$data$tmin + wacsdata$data$trange
      wacsdata$data = cbind(wacsdata$data,tmax)
      tmax = wacssimul$sim$tmin + wacssimul$sim$trange
      wacssimul$sim = cbind(wacssimul$sim,tmax)
    }
  }
  
  #
  # Different types of validations
  #
  
  #############################   VALIDATING RAIN
  if (what == "Rain") {
    if (is.null(wacsdata) || is.null(wacspar)) {
      stop ("[WACSvalid] for 'Rain' you must provide wacsdata and wacspar")
    }  
    res = wacsvalid.Rain(wacsdata,wacspar)
  }
  
  #############################   VALIDATING ONE SIMULATION RUN
  if (what == "Sim"){
    if (is.null(wacsdata) || is.null(wacssimul) ||  is.null(wacspar) || is.null(varname)) {
      stop ("[WACSvalid] for 'Sim' you must provide wacsdata, wacssimul, wacspar
               and varname");
    }
    res = wacsvalid.Sim(wacsdata,wacspar,wacssimul,varname)
  }
  
  #############################   VALIDATING Means and Standard Deviations
  if (what == "MeanSd") {
    if (is.null(wacsdata) || is.null(wacssimul) || is.null(varname)) {
      stop ("[WACSvalid] for 'MeanSd' you must provide wacsdata, wacssimul
               and varname");
    }
    res = wacsvalid.MeanSd(wacsdata,wacssimul,varname)
  }
  
  #############################   VALIDATING Bivariate correlation
  if (what == "BiVar") {
    if (is.null(wacsdata) || is.null(wacssimul) || is.null(varname) || is.null(varname2))  {
      stop ("[WACSvalid] for 'BiVar' you must provide wacsdata, 
            wacssimul, varname and varname2")}
    if ((varname2=="tmax") && wacsdata$Trange){
      tmax = wacsdata$data$tmin + wacsdata$data$trange
      wacsdata$data = cbind(wacsdata$data,tmax)
      tmax = wacssimul$sim$tmin + wacssimul$sim$trange
      wacssimul$sim = cbind(wacssimul$sim,tmax)
    } 
    res=wacsvalid.BiVar(wacsdata,wacssimul,varname,varname2)
  }
  
  #############################   VALIDATING Temporal Correlations
  if (what == "CorTemp") {
      if (is.null(wacsdata) || is.null(wacssimul) || is.null(varname)) {
          stop ("[WACSvalid] for 'CorTemp' you must provide wacsdata, 
                  wacssimul and varname");
      }
      res = wacsvalid.CorTemp(wacsdata,wacssimul,varname)
  }    
  
  #############################   VALIDATING Sum above base
  if (what == "SumBase") {
    if ( is.null(wacsdata) || is.null(wacssimul) || is.null(varname) || is.null(base) || is.null(months) ) {
      stop ("[WACSvalid] for 'SumBase' you must provide wacsdata, 
                  wacssimul, varname, base, months");
    }
    res = wacsvalid.SumBase(wacsdata,wacssimul,varname,base,months)
  }   
  
  #############################   VALIDATING Persistance
  if (what == "Persistence") {
    if ( is.null(wacsdata) || is.null(wacssimul) || is.null(varname) || is.null(base) || is.null(months) || is.null(above)) {
       stop ("[WACSvalid] for 'Persistence' you must provide wacsdata,
                        wacssimul, varname, base, above, months");
    }
    res = wacsvalid.Persistence(wacsdata,wacssimul,varname,base,above,months)
  }
  res$labels = c("Observed","Simulated")
  return(res)
}
