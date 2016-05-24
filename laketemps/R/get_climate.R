#'@title get climate data according to lake name from GLTC dataset
#'@description Climate data compiled for the Global Lake Temperature Collaboration can be 
#'accessed and returned as data.frames. See associated publication and references therein for 
#'details including units and data provenance. 
#'@param lake_name a valid name of a lake in the GLTC dataset (see \code{\link{get_lake_names}}). 
#'\code{lake_name} is case insensitive. 
#'@param types name for the climate data (see \code{\link{get_climate_names}}). 
#'\code{types} is case insensitive. 
#'All data are returned if this argument is missing.
#'@return a lake data.frame, or an empty data.frame if no data exist
#'@export
#'@seealso \code{\link{get_climate_names}}, \code{\link{get_lake_names}}, \code{\link{get_surface_temps}}
#'@references
#'Sharma, Sapna; Gray, Derek; Read, Jordan; O'Reilly, Catherine; Schneider, Philipp; Qudrat, Anam; Gries, Corinna;
#'Stefanoff, Samantha; Hampton, Stephanie; Hook, Simon; Lenters, John; Livingstone, David; McIntyre, Peter; Adrian, Rita;
#'Allan, Mathew; Anneville, Orlane; Arvola, Lauri; Austin, Jay; Bailey, John; Baron, Jill; Brookes, Justin; Chen, Yuwei;
#'Daly, Robert; Dokulil, Martin; Dong, Bo; Ewing, Kye; de Eyto, Elvira; Hamilton, David; Havens, Karl; Haydon, Shane;
#'Hetzenauer, Harald; Heneberry, Joceylene; Hetherington, Amy; Higgins, Scott; Hixcon, Eric; Izmest'eva, Lyubov; Jones,
#'Benjamin; Kangur, Kullli; Kasprzak, Peter; Koster, Olivier; Kraemer, Benjamin; Kumagai, Michio; Kuusisto, Esko;
#'Leshkevich, George; May, Linda; MacIntyre, Sally; Mueller-Navarra, Doerthe; Naumenko, Mikhail; Noges, Peeter; Noges,
#'Tiina; Niederhauser, Pius; North, Ryan; Paterson, Andrew; Plisnier, Pierre-Denis; Rigosi, Anna; Rimmer, Alon; Rogora,
#'Michela; Rudstram, Lars; Rusak, James; Salmaso, Nico; Samal, Nihar; Schindler, Daniel; Schladow, Geoffrey; Schmidt, Silke;
#'Schultz, Tracey; Silow, Eugene; Straile, Dietmar; Teubner, Katrin; Verburg, Piet; Voutilainen, Ari; Watkinson, Andrew;
#'Weyhenmeyer, Gesa; Williamson, Craig; Woo, Kara (2014): Globally distributed lake surface water temperatures collected in
#'situ and by satellites; 1985-2009. Long Term Ecological Research Network.
#'http://dx.doi.org/10.6073/pasta/379a6cebee50119df2575c469aba19c5
#'@examples
#'get_climate_names()
#'get_climate('Victoria', types = c('Radiation.Shortwave.Summer', 'Air.Temp.Mean.Summer.NCEP'))
#'get_climate('Mendota')
#'get_climate('mendota','radiation.shortwave.summer')
get_climate <- function(lake_name, types){
  

  
  if (missing(types)){
    types <- get_climate_names()
  }
  check_climate(types)
  
  climate_data <- subset_lake_data(lake_name, types)
  return(climate_data)
  
}

