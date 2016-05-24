#' biolflor_lookup
#' 
#' This dataframe is used to retrieve species URLs from the BiolFlor website
#' (\samp{http://www2.ufz.de/biolflor/index.jsp}). 
#' 
#' @name biolflor_lookup
#' @docType data
#' @format A data frame with 3688 observations on the following 8 variables.
#' \describe{
#' \item{submittedname}{a character vector}
#' \item{acceptedname}{a character vector}
#' \item{sourceid}{a character vector}
#' \item{score}{a character vector}
#' \item{matchedname}{a character vector}
#' \item{annotations}{a character vector}
#' \item{V1}{a character vector}
#' \item{V2}{a character vector} }
#' @references Please use the follow citation ay time you use data derived from Biolflor:
#' BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland.
#' Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)
#' @keywords datasets
#' @examples \dontrun{
#' data(biolflor_lookup)
#' }
#' @usage  data(biolflor_lookup)
NULL

#' biolflor_check
#' 
#' A lookup dataframe for checking species names on the Biolflor website
#' confronting them with the \code{tnrs} function
#' 
#' 
#' @name biolflor_check
#' @docType data
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @format A data frame with 3688 observations on the following 8 variables.
#' \describe{
#' \item{submittedname}{name of plant species submitted to \code{tnrs}}
#' \item{acceptedname}{accepted name according to \code{tnrs}}
#' \item{sourceid}{sourceid according to \code{tnrs}}
#' \item{score}{score of matching species names according to \code{tnrs}}
#' \item{matchedname}{matched name according to \code{tnrs}}
#' \item{annotations}{plat species authors}
#' \item{uri}{url for the plant species on the \code{tnrs} website}
#'  }
#' @references Please use the follow citation ay time you use data derived from Biolflor:
#' BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland.
#' Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)
#' @keywords dataframe
#' @examples \dontrun{
#' #data(biolflor_check)
#' }
#' 
NULL

#' list_of_traits_Biolflor
#'
#' a vector containing traits that can be downloaded from Biolflor
#'
#'
#' @name list_of_traits_Biolflor
#' @docType data
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @format A vector of plant traits
#' @references Please use the follow citation ay time you use data derived from Biolflor:
#' BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland.
#' Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)
#' @keywords list_of_traits_Biolflor
#' @examples \dontrun{
#' #data(list_of_traits_Biolflor)
#' }
NULL

#' traits_special_Biolflor
#' 
#' a vector containing traits that can be downloaded from Biolflor for which some special Xpaht rules must be applied
#'
#' @name traits_special_Biolflor
#' @docType data
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @format A vector of plant traits
#' @references Please use the follow citation ay time you use data derived from Biolflor:
#' BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland.
#' Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)
#' @keywords traits_special_Biolflor
#' @examples \dontrun{
#' #data(traits_special_Biolflor)
#' }
NULL

#' traits_pollen_Biolflor
#'
#' a vector containing traits that can be downloaded from Biolflor
#'
#'
#' @name traits_pollen_Biolflor
#' @docType data
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @format A vector of plant traits
#' @references Please use the follow citation ay time you use data derived from Biolflor:
#' BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland.
#' Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)
#' @keywords biolflor
#' @examples \dontrun{
#' #data(traits_pollen_Biolflor)
#' }
NULL


#' \bold{ECOFLORA_df}: local lookup table for Ecoflora species' url(s)
#' 
#' This dataset is not ment to be directly accessed by the final user. It is
#' imported by the \code{ecoflora()} function to extrapolate the correct
#' \code{URL} for each species of interest and download the corresponding
#' functional traits.
#' This dataset is used as a lookup table from the \code{ecoflora()} function.
#' 
#' @name ECOFLORA_df
#' @docType data
#' @format A data frame containing URL for the 1879 species contained in the
#' Ecoflora web database.  \describe{
#' \item{species}{a vector containing the species names as defined on Ecoflora website}
#' \item{web_link}{a character vector containing the URL of each species trait web page}
#' \item{acceptedname}{a character vector containing the accepted name according to TNRS}
#' \item{sourceid}{a character vector containing the source used by the \code{taxize::tnrs} function}
#' \item{score}{a numeric vector containing the score obtained by \code{taxize::tnrs} function}
#' \item{matchedname}{a character vector containing the matched names by \code{taxize::tnrs}}
#' \item{annotations}{a character vector}
#' \item{uri}{a character vector containing Ecoflora-URL for each species}}
#' @references Please alwasy cite the database according to:
#' 
#' Fitter, A . H. and Peat , H. J., 1994, The Ecological Flora Database, J.
#' Ecol., 82, 415-425.
#' 
#' @keywords datasets
#' @examples \dontrun{
#' data(ECOFLORA_df)
#' }
NULL

#' Set of functional traits to be retrieved by Ecoflora.
#' 
#' \code{traits_eco} defines a list containg pairs in the form
#' \emph{short_name_of_the_trait} = \emph{corresponding_code_in_Ecoflora }
#' At the moment the package does not download all the traits available at
#' Ecoflora; curious users can expand the number of downloadable traits
#' simply extending the list with other 'pairs' (take care of using the right
#' \code{Ecoflora codes} as described in
#' \samp{http://www.ecoflora.co.uk/search_plantchar.php}).
#' 
#' @name traits_eco
#' @docType data
#' @format The format is: a list of the following 17 elements, where
#' each element is a pair of the form "traits":"code used in Ecoflora HTML code":\describe{
#'  \item{height_max }{: num 3.05}
#'  \item{height_min }{: num 3.06}
#'  \item{leaf_area }{: num 3.17}
#'  \item{leaf_longevity }{: num 3.22}
#'  \item{Photosynthetic_pathway }{: num 4.02}
#'  \item{life_form }{: num 5.01}
#'  \item{Vegetative_reprod_method}{: num 5.05}
#'  \item{Flowering_earliest_month}{: num 5.07}
#'  \item{Flowering_latest_month }{: num 5.08}
#'  \item{Pollen_vector }{: num 5.15}
#'  \item{Seed_weight_mean }{: num 5.34}
#'  \item{Method_of_propagation }{: num 5.52}
#'  \item{Ellenberg_light_Eco }{: num 7.14}
#'  \item{Ellenberg_moisture_Eco }{: num 7.15}
#'  \item{Ellenberg_pH_Eco }{: num 7.16}
#'  \item{Ellenberg_nitrogen_Eco }{: num 7.17}
#'  \item{Ellenberg_salt_Eco }{: num 7.18}
#' }
#' @docType data
#' @references Fitter, A . H. and Peat , H. J., 1994, The Ecological Flora Database,
#' J. Ecol., 82, 415-425.  \samp{http://www.ecoflora.co.uk}
#' @keywords datasets
#' @examples \dontrun{
#' data(traits_eco)
#' }
NULL


###* LEDA

#' \bold{LEDA_df} : local (partial) copy of the LEDA traitbase website.
#' 
#' This data.frame includes a local copy of (some) of the parameters available
#' on the LEDA website; up to now only few functional traits were taken into
#' account.  Hopefully the database will be extendend in the following releases
#' of the package.  Please refer to
#' \samp{http://www.leda-traitbase.org/LEDAportal/plantTraits.jsp} for a
#' detailed explanation of the traits availabe (definitions, units of measure,
#' etc...) in LEDA.
#' 
#' The original datasets are available as .TXT files in they raw form at
#' \samp{http://www.leda-traitbase.org/LEDAportal/data_files.jsp}; what is found
#' here is a selection of a few traits.  A higher number of traits from LEDA
#' may be included in the following releases of the package.
#' 
#' 
#' @name LEDA_df
#' @docType data
#' @format A data frame containing traits data for 8309 species.  \describe{
#' \item{SBS.name}{species name as used in LEDA}
#' \item{canopy_height.m}{: canopy hight expressed in meters}
#' \item{mean.SLA..mm.2.mg.}{: mean value for specific leaf area expressed in \emph{mm^2 mm^-1}}
#' \item{mean.SM..mg.}{: seed mass in \code{g}}
#' \item{plant.growth.form}{: plant growth form }
#' \item{dispersal.type}{: dispersal type categories}
#' \item{acceptedname}{: accepted name according to \emph{TNRS}}
#' \item{seed_longevity}{: longevity of the seedbank}
#' \item{sourceid}{: source of species name used by \code{taxize::tnrs}}
#' \item{score}{: score obtained by  \code{taxize::tnrs}}
#' \item{matchedname}{: matched name obtained by \code{taxize::tnrs}}
#' \item{annotations}{: }
#' \item{uri}{: reference \code{url} for the species} }
#' @references Please cite the following reference any time you use data retrieved from the
#' LEDA traitbase (citation reported at
#' \samp{http://www.leda-traitbase.org/LEDAportal/citation.jsp}):
#' 
#' Kleyer, M., Bekker, R.M., Knevel, I.C., Bakker, J.P, Thompson, K.,
#' Sonnenschein, M., Poschlod, P., Van Groenendael, J.M., Klimes, L.,
#' Klimesova, J., Klotz, S., Rusch, G.M., Hermy, M., Adriaens, D., Boedeltje,
#' G., Bossuyt, B., Dannemann, A., Endels, P., Götzenberger, L., Hodgson, J.G.,
#' Jackel, A-K., Kühn, I., Kunzmann, D., Ozinga, W.A., Römermann, C., Stadler,
#' M., Schlegelmilch, J., Steendam, H.J., Tackenberg, O., Wilmann, B.,
#' Cornelissen, J.H.C., Eriksson, O., Garnier, E., Peco, B. (2008): The LEDA
#' Traitbase: A database of life-history traits of Northwest European flora.
#' Journal of Ecology 96: 1266-1274.
#' @keywords datasets
#' @examples \dontrun{
#' data(LEDA_df)
#' }
#' 
NULL


#' TR8: a tool for retrieving functional traits data for plant species.
#' 
#' This package provide a set of functions for extracting traits data for plant
#' species from various sources.
#' 
#' \tabular{ll}{ Package: \tab TR8\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2014-02-27\cr License: \tab GPL>=2\cr Depends: \tab XML,
#' RCurl, plyr, taxize\cr }
#' 
#' The easiest way of using the package is through the \code{tr8()} function,
#' which accepts a vector of plant species names and returns a data frame
#' containg traits data which have been found in the various sources.
#' 
#' @name TR8-package
#' @aliases TR8-package TR8
#' @docType package
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @keywords package, functional traits
#' @import XML RCurl plyr methods
#' @examples \dontrun{
#' ## store the resulting data.frame in \code{My_traits}
#' 
#' #My_traits<-tr8(c("Abutilon theophrasti","Trifolium campestre"))
#' }
NULL

#' AMF inoculation potential
#'
#'
#' 
#' @format a dataframe of 240 rws with the following 2 variables \describe{
#'
#'   \item{species}{name of plant species}
#'   \item{support.AMF}{score which assess the "potentiality for the species to host AMF"}
#' 
#' }
#' @name AMF_giovannetti 
#' @docType data
#' @references Giovannetti E., Avio L., Sbrana C. (personal communication)
#' @seealso \code{help(retrieve_amf)}
#' @keywords datasets
#' @examples \dontrun{
#'#' data(AMF_giovannetti)
#' }
#'
NULL


#' pignatti
#' 
#' \code{pignatti} is a dataframe containing traits data
#' for italian species
#' 
#' @name pignatti
#' @docType data
#' @format The format is: a list of the following 17 elements, where
#' each element is a pair of the form "traits":"code used in Ecoflora HTML code":\describe{
#'  \item{Specie.Pignatti }{: species name in the original dataset }
#'  \item{numero          }{: numeric code }
#'  \item{codice          }{: numeric code }
#'  \item{nome.scientifico}{: scientific name with authors }
#'  \item{forma_biologica }{: life form}
#'  \item{corotipo        }{: distribution of species }
#'  \item{L               }{: Ellenberg value for light}
#'  \item{T               }{: Ellenberg value for temperature}
#'  \item{C               }{: Ellenberg value for continentality}
#'  \item{U               }{: Ellenberg value for soil humidity }
#'  \item{R               }{: Ellenberg value for soil pH}
#'  \item{N               }{: Ellenberg value for Nutrients in the soil }
#'  \item{S               }{: Ellenberg value for soil salinity}
#'  \item{Name.tnrs       }{: species name according to tnrs }
#' }
#' @docType data
#' @references Pignatti S., Menegoni P., Pietrosanti S., 2005, Biondicazione attraverso le piante vascolari. Valori di indicazione secondo Ellenberg (Zeigerwerte) per le specie della Flora d'Italia. Braun-Blanquetia 39, Camerino, pp.  97.
#' @keywords datasets
#' @examples \dontrun{
#' data(pignatti)
#' }
NULL


#' myco
#' 
#' A lookup dataframe for checking species names on the Biolflor website
#' confronting them with the \code{tnrs} function
#' 
#' 
#' @name myco
#' @docType data
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @format A data frame with 2970 observations on the following 2 variables.
#' \describe{
#' \item{species}{name of plant species}
#' \item{Myco_infection}{score of AMF infection}
#'  }
#' @references  Asem A. Akhmetzhanova, Nadejda A. Soudzilovskaia, Vladimir G. Onipchenko,
##' Will K. Cornwell, Vladimir A. Agafonov, Ivan A. Selivanov, and Johannes H. C. Cornelissen.
##' 2012. A rediscovered treasure: mycorrhizal intensity database for 3000 vascular plants
##' species across the former Soviet Union. Ecology 93:689.
##' \samp{http://esapubs.org/Archive/ecol/E093/059/default.htm}
##'
##'
#' @keywords dataframe
#' @examples \dontrun{
#' #data(myco)
#' }
#' 
NULL


#' column_list
#' 
#' A list containing a brief description of traits data retrieved
#' by the various databases
#' 
#' 
#' @name column_list
#' @docType data
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @format A list of 40 elements: \describe{
#' \item{height_max}{: c(h_max,maximum height,Ecoflora)}
#' \item{height_min}{: c(h_min,minimum height,Ecoflora)}
#' \item{leaf_area}{: c(le_area,leaf area,Ecoflora)}
#' \item{leaf_longevity}{: c(le_long,leaf longevity,Ecoflora)}
#' \item{Photosynthetic_pathway}{: c(phot_path,photosynthetic pathway,Ecoflora)}
#' \item{life_form}{: c(li_form,life form,Ecoflora)}
#' \item{Vegetative_reprod_method}{: c(reprod_meth,vegetative reprod method,Ecoflora)}
#' \item{Flowering_earliest_month}{: c(flw_early,flowering earliest month,Ecoflora)}
#' \item{Flowering_latest_month}{: c(flw_late,Flowering_latest_month,Ecoflora)}
#' \item{Pollen_vector}{: c(poll_vect,pollen vector,Ecoflora)}
#' \item{Seed_weight_mean}{: c(seed_wght,seed weight mean,Ecoflora)}
#' \item{Method_of_propagation}{: c(propag,Method of propagation,Ecoflora)}
#' \item{Ellenberg_light_Eco}{: c(ell_light_uk,Ellenberg light,Ecoflora)}
#' \item{Ellenberg_moisture_Eco}{: c(ell_moist_uk,Ellenberg moisture,Ecoflora)}
#' \item{Ellenberg_pH_Eco}{: c(ell_pH_uk,Ellenberg pH,Ecoflora)}
#' \item{Ellenberg_nitrogen_Eco}{: c(ell_N,Ellenberg nitrogen,Ecoflora)}
#' \item{Ellenberg_salt_Eco}{: c(ell_S,Ellenberg salt,Ecoflora)}
#' \item{canopy_height.m}{: c(c_height,canopy height in [m],LEDA)}
#' \item{mean.SLA..mm.2.mg.}{: c(SLA,mean specific leaf area [mm^2 m-1],LEDA)}
#' \item{mean.SM..mg.}{: c(SM,mean seed mass,LEDA)}
#' \item{plant.growth.form}{: c(g_form,plant growth form,LEDA)}
#' \item{dispersal.type}{: c(disp_type,plant dispersal type,LEDA)}
#' \item{Life.form}{: c(li_form_B,life form,BiolFlor)}
#' \item{Life.span}{: c(li_span,life span,BiolFlor)}
#' \item{Rosettes}{: c(ros,Rosettes,BiolFlor)}
#' \item{Type.of.reproduction}{: c(reprod_B,Type of reproduction,BiolFlor)}
#' \item{Strategy.type}{: c(strategy,Strategy type,BiolFlor)}
#' \item{Pollen.vector}{: c(poll_vect_B,Pollen vector,BiolFlor)}
#' \item{Flower.class.after.MUELLER}{: c(flw_muell,Flower class after MUELLER,BiolFlor)}
#' \item{Begin.of.flowering..months.}{: c(flw_beg,Begin of flowering months.,BiolFlor)}
#' \item{End.of.flowering..months.}{: c(flw_end,End of flowering months,BiolFlor)}
#' \item{Duration.of.flowering..months.}{: c(flw_dur,Duration of flowering months,BiolFlor)}
#' \item{Number.of.flowering.phases}{: c(flw_ph,Number of flowering phases,BiolFlor)}
#' \item{L}{: c(ell_L_it,ellenberg value for light in Italy,Pignatti et al.)}
#' \item{T}{: c(ell_T_it,ellenberg value for temperature in Italy,Pignatti et al.)}
#' \item{C}{: c(ell_C_it,ellenberg value for continentality in Italy,Pignatti et al.)}
#' \item{U}{: c(ell_U_it,ellenberg value for humidity in Italy,Pignatti et al.)}
#' \item{R}{: c(ell_R_it,ellenberg value for soil reaction in Italy,Pignatti et al.)}
#' \item{N}{: c(ell_N_it,ellenberg value for nitrogen in Italy,Pignatti et al.)}
#' \item{S}{: c(ell_S_it,ellenberg value for salinity in Italy,Pignatti et al.)}
#'  }
#' @references  Asem A. Akhmetzhanova, Nadejda A. Soudzilovskaia, Vladimir G. Onipchenko,
##' Will K. Cornwell, Vladimir A. Agafonov, Ivan A. Selivanov, and Johannes H. C. Cornelissen.
##' 2012. A rediscovered treasure: mycorrhizal intensity database for 3000 vascular plants
##' species across the former Soviet Union. Ecology 93:689.
##' \samp{http://esapubs.org/Archive/ecol/E093/059/default.htm}
##'
##'
#' @keywords dataframe
#' @examples \dontrun{
#' #data(myco)
#' }
#' 
NULL

