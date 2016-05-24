#' Simple Event Table
#'
#' A simple, hypothetical event table.
#'
#' \itemize{
#'   \item from, to - endpoint positions
#'   \item x, y, z - numeric variables
#'   \item factor - a factor variable
#' }
#'
#' @format A data frame with 11 rows and 7 variables.
#' @name simple
NULL

#' Elwha River Survey
#'
#' An event table containing the results of a survey of the Elwha River (Washington, USA) carried out in August-September 2008. Both physical variables and fish counts were collected.
#'
#' \itemize{
#'   \item from, to - distance upstream from the river mouth [km]
#'   \item unit.length - unit length [m]
#'   \item unit.type - unit type (P = pool, GP = glide-like pool, G = glide, GR = glide-like riffle, R = riffle)
#'   \item channel.type - channel type (1 = main, 2 = secondary)
#'   \item mean.depth - mean depth [m]
#'   \item max.depth - max depth [m]
#'   \item mean.width - mean wetted width [m]
#'   \item bedrock - bedrock substrate [\%]
#'   \item boulder - boulder substrate [\%]
#'   \item cobble - cobble substrate [\%]
#'   \item gravel - gravel substrate [\%]
#'   \item sand - sand substrate [\%]
#'   \item silt - silt substrate [\%]
#'   \item overhang.cover - channel banks with overhanging vegetation [\%]
#'   \item boulder.cover - channel area covered by boulders [\%]
#'   \item jams - number of log jams
#'   \item jam.area - total area of log jams [m^2]
#'   \item SACO.10/20/30/40/total - number of Bull Trout (\emph{Salvelinus confluentus}) sized 10 - 20 cm / 20 - 30 cm / 30 - 40 cm / > 40 cm / total, respectively.
#'   \item ONXX.10/20/30/40/total - number of trout (\emph{Oncorhynchus sp.}) sized 10 - 20 cm / 20 - 30 cm / 30 - 40 cm / > 40 cm / total, respectively.
#'   \item SAFO - number of Brook Trout (\emph{Salvelinus fontinalis})
#'   \item ONTS - number of Chinook Salmon (\emph{Oncorhynchus tshawytscha})
#'   \item ONNE - number of Sockeye Salmon (\emph{Oncorhynchus nerka})
#'   \item LATR - number of Pacific Lamprey (\emph{Lampetra tridentata})
#'   \item ONKI - number of Coho Salmon (\emph{Oncorhynchus kisutch})
#' }
#'
#' @format A data frame with 249 rows and 33 variables.
#' @source Brenkman, S. J., J. J. Duda, C. E. Torgersen, E. Z. Welty, G. R. Pess, R. Peters, and M. L. McHenry. 2012. A riverscape perspective of Pacific salmonids and aquatic habitats prior to large-scale dam removal in the Elwha River, Washington, USA. Fisheries Management and Ecology 19:36-53. DOI: \href{http://dx.doi.org/10.1111/j.1365-2400.2011.00815.x}{10.1111/j.1365-2400.2011.00815.x}
#' @name elwha
NULL

#' Quinault River Survey
#'
#' An event table containing the results of a survey of the Quinault River (Washington, USA) in August 2009. Both physical variables and fish counts were collected.
#'
#' \itemize{
#'   \item from, to - distance upstream from the river mouth [km]
#'   \item altitude - mean elevation above sea level [m]
#'   \item channel.type - channel type (1 = main, 2 = secondary)
#'   \item unit.type - unit type (P = pool, GP = glide-like pool, GR = glide-like riffle, R = riffle)
#'   \item unit.length - unit length [m]
#'   \item mean.width - mean wetted width [m]
#'   \item mean.depth - mean depth [m]
#'   \item max.depth - max depth [m]
#'   \item overhang.cover - channel banks with overhanging vegetation [\%]
#'   \item boulder.cover - channel area covered by boulders [\%]
#'   \item jams - number of log jams
#'   \item jam.area - total area of log jams [m^2]
#'   \item SACO.10/20/30/50/total - number of Bull Trout (\emph{Salvelinus confluentus}) sized 10 - 20 cm / 20 - 30 cm / 30 - 50 cm / > 50 cm / total, respectively.
#'   \item ONXX.10/20/30/total - number of trout (\emph{Oncorhynchus sp.}) sized 10 - 20 cm / 20 - 30 cm / > 30 cm / total, respectively.
#'   \item PRWI - number of Mountain Whitefish (\emph{Prosopium williamsoni})
#'   \item ONTS - number of Chinook Salmon (\emph{Oncorhynchus tshawytscha})
#'   \item ONMY - number of Rainbow Trout (\emph{Oncorhynchus mykiss})
#'   \item ONKI - number of Coho Salmon (\emph{Oncorhynchus kisutch})
#'   \item ONNE - number of Sockeye Salmon (\emph{Oncorhynchus nerka})
#'   \item ONGO - number of Pink Salmon (\emph{Oncorhynchus gorbushcha})
#'   \item ONKE - number of Chum Salmon (\emph{Oncorhynchus keta})
#'   \item CAMA - number of Largescale Sucker (\emph{Catostomus macrocheilus})
#'   \item LATR - number of Pacific Lamprey (\emph{Lampetra tridentata})
#' }
#'
#' @format A data frame with 363 rows and 31 variables.
#' @source Samuel J. Brenkman (National Park Service, Olympic National Park, Washington, USA), unpublished data.
#' @name quinault
NULL

#' Dungeness River (NetMap)
#' 
#' NetMap (\href{http://www.terrainworks.com/}{terrainworks.com}) output for the entire fluvial network of the Dungeness River (Washington, USA). NetMap employs digital elevation models to generate detailed river networks and compute biophysical variables for spatially continuous hydrologic units throughout the networks.
#' 
#' \itemize{
#'   \item CHAN_ID - channel identifier (1 = mainstem, all others are tributaries)
#'   \item OUT_DIST - distance upstream from the river mouth to the downstream end of the unit [km]
#'   \item LENGTH_M - unit length [m]
#'   \item ... - see NetMap's \href{http://www.netmaptools.org/Pages/NetMapHelp/master_attribute_list.htm}{Master Attribute List}
#' }
#'
#' @format A data frame with 16,616 rows and 47 variables.
#' @source \url{http://www.terrainworks.com/netmap-demo-tools-download}
#' @name netmap
NULL

#' Fish Movements
#' 
#' A pair of event tables (in a list) documenting the movements of tagged Coho Salmon (\emph{Oncorhynchus kisutch}) in Bear Creek (Southwest Alaska, USA) for 29 July - 19 August 2008. Table \code{motion} lists individual fish residence time intervals in each of three stream regions, while table \code{origin} lists the study-wide residence time of each fish and the stream region in which the fish was first tagged.
#' 
#' \itemize{
#'   \item from, to - start and end times as seconds since 1970-01-01 UTC (POSIXct)
#'   \item fish.id - unique identifer for each fish
#'   \item region - stream region (1 = 0 - 930 m, a cold downstream region with abundant and spawning sockeye salmon; 2 = 930 - 1360 m, a cold middle region with few if any sockeye salmon; 3 = > 1360 m, a warm upstream region where sockeye salmon were absent)
#' }
#'
#' @format Two data frames \code{motion} and \code{origin} with 1,140 rows and 149 rows of 4 variables, respectively.
#' @source Armstrong, J. B., D. E. Schindler, C. P. Ruff, G. T. Brooks, K. E. Bentley, and C. E. Torgersen. 2013. Diel horizontal migration in streams: juvenile fish exploit spatial heterogeneity in thermal and trophic resources. Ecology 94:2066-2075. DOI: \href{http://dx.doi.org/10.1890/12-1200.1}{10.1890/12-1200.1}
#' @name fishmotion
NULL