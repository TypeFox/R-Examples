#'@title get summer surface lake temperature data from GLTC dataset
#'@description get summer lake surface temperatures for the Global Lake Temperature Collaboration dataset. 
#'All temperatures are in degrees C. 
#'@param lake_name a valid name of a lake in the GLTC dataset (see \code{\link{get_lake_names}}). 
#'\code{lake_name} is case insensitive. 
#'@param type source for the data. Either "Lake.Temp.Summer.InSitu" or "Lake.Temp.Summer.Satellite". 
#'if missing, both sources will be used.
#'@return a lake data.frame (empty data.frame) if no data exist. Temperatures are returned in degrees C.
#'@seealso \code{\link{get_lake_names}}, \code{\link{get_climate}}
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
#'get_surface_temps('Victoria','Lake.Temp.Summer.Satellite')
#'get_surface_temps('Mendota','Lake.Temp.Summer.InSitu')
#'
#'# - expect no satellite data for Lake Mendota:
#'get_surface_temps('Mendota','Lake.Temp.Summer.Satellite')
#'
#'# retrieve from a lake site with multiple surface temperature sources:
#'get_surface_temps('Tahoe.Mid.Lake')
#'# is the same as:
#'get_surface_temps('Tahoe.Mid.Lake',c('Lake.Temp.Summer.Satellite', 'Lake.Temp.Summer.InSitu'))
#'@export
get_surface_temps <- function(lake_name, type){
  
  both = temp_types()
  
  if (missing(type)){
    type = both
  }
  
  if (!all(type %in% both)){stop(paste0('type=', type, ' not recognized'))}
  
  surf_temps <- subset_lake_data(lake_name, types = type)
  return(surf_temps)
  
}

