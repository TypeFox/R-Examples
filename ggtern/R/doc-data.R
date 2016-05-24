#' @title Elkin and Groves Feldspar Data
#' 
#' @description Data relating to Elkins and Groves Feldspar Data, the following datasets include 
#' the experimental data and sample raster data from one of the images in the 
#' referenced paper.
#' \code{Feldspar} - Experimental Data
#' \code{FeldsparRaster} - Raster Data for Fig. 6.
#' 
#' @references 
#' Elkins, L. T. & Grove, T. L. 
#' Ternary Feldspar Experiments and Thermodynamic Models 
#' American Mineralogist, Mineral Soc America, 1990, 75, 544-559
#' @docType data
#' @usage 
#' #Experimental Data
#' data(Feldspar)
#' 
#' #Raster data
#' data(FeldsparRaster)
#' @format 
#' \code{Feldsdpar} - One (1) row per Feldspar composition, \code{FeldsdparRaster} - Raster Matrix
#' @examples
#' #Plot Felspar Data
#' data(Feldspar)
#' summary(Feldspar)
#' ggtern(data=Feldspar,aes(x=An,y=Ab,z=Or)) + geom_point()
#' 
#' # Plot Feldspar data and Underlying Raster Image
#' data(FeldsparRaster)
#' ggtern(Feldspar,aes(Ab,An,Or)) + 
#' theme_rgbw() + 
#' annotation_raster_tern(FeldsparRaster,xmin=0,xmax=1,ymin=0,ymax=1) +
#' geom_point(size=5,aes(shape=Feldspar,fill=Feldspar),color='black') +
#' scale_shape_manual(values=c(21,24)) +
#' labs(title="Demonstration of Raster Annotation")
#' @seealso \link[=data]{Data}
#' @name data_sets_Feldspar
#' @rdname data_sets_Feldspar
#' @aliases Feldspar FeldsparRaster
#' @author Nicholas Hamilton
NULL

#' @title USDA Textural Classification Data
#' 
#' @description This dataset was issued by the United States Department of Agriculture (USDA) 
#' in the form of a ternary diagram, this original ternary diagram has been converted to numerical data 
#' and included here.
#' @docType data
#' @usage data(USDA)
#' @format 1row per point, many points per classification representing the extremes of the area.
#' @source Soil Mechanics Level 1, Module 3, USDA Textural Classification Study Guide
#' @author United States Department of Agriculture (USDA)
#' @seealso \link[=data]{ggtern datasets}
#' @examples
#' #Load the Libraries
#' library(ggtern)
#' library(plyr)
#' #Load the Data.
#' data(USDA)
#' #Put tile labels at the midpoint of each tile.
#' USDA.LAB <- ddply(USDA,"Label",function(df){apply(df[,1:3],2,mean)})
#' #Tweak
#' USDA.LAB$Angle=0; USDA.LAB$Angle[which(USDA.LAB$Label == "Loamy Sand")] = -35
#' #Construct the plot.
#' ggtern(data=USDA,aes(Sand,Clay,Silt,color=Label,fill=Label)) +
#' geom_polygon(alpha=0.75,size=0.5,color="black") +
#' geom_mask() +  
#' geom_text(data=USDA.LAB,aes(label=Label,angle=Angle),color="black",size=3.5) +
#' theme_rgbw() + 
#' theme_showsecondary() +
#' theme_showarrows() +
#' weight_percent() + guides(fill='none') + 
#' theme_legend_position("topleft")
#' labs(title="USDA Textural Classification Chart",fill="Textural Class",color="Textural Class")
#' @name data_sets_USDA
#' @rdname data_sets_USDA
#' @aliases USDA
#' @author Nicholas Hamilton
NULL



#' Grantham and Valbel Rock Fragment Data
#' 
#' \strong{ABSTRACT:} Chemical weathering influences the detrital composition of sand-size sediment derived from source 
#' areas subject to different amounts of precipitation in the Coweeta Basin, North Carolina. Of the grain types 
#' studied, rock fragments are most sensitive to chemical degradation; therefore, their abundance is the best 
#' indicator of cumulative weathering effects. Destruction of sand-size rock fragments by chemical weathering 
#' is a function of both the intensity and duration of chemical weathering experienced by grains in regoliths 
#' of the source area. In the Coweeta Basin, the intensity of chemical weathering is directly related to the 
#' climate via effective precipitation in individual subbasins, whereas the duration of chemical weathering is 
#' inversely related to the relief ratio of the watershe . Therefore, soils in watersheds with low-relief 
#' ratios and high discharge per unit area experience the most extensive chemical weathering, and sediments 
#' derived from these watersheds contain the lowest percentage of rock fragments. The effects of climate alone 
#' cannot explain the systematic variation of rock fragment abundance in sediments from the Coweeta Basin. 
#' The compositional imprint left on these sediments by chemical weathering is a function of both climate and 
#' topographic slope in the sediment source area.
#' @docType data
#' @references Grantham, Jeremy Hummon, and Michael Anthony Velbel. 
#' "The influence of climate and topography on rock-fragment abundance in modern fluvial sands of the southern 
#' Blue Ridge Mountains, North Carolina." Journal of Sedimentary Research 58.2 (1988).
#' @usage data(Fragments)
#' @format 1row per point, Each point contains data on the following:
#' \enumerate{ 
#' \item \strong{Watershed}: By id: 2, 10, 34, 41, 13, 27, 32 or 37, 
#' \item \strong{Position}: By name: Tallulah or Coweeta, 
#' \item \strong{CCWI}: The Cumulative Chemical Weathering Index: numeric
#' \item \strong{Precipitation}: Average Annual Precipitation, numeric 
#' \item \strong{Discharge}: Annual Average Discharge, numeric 
#' \item \strong{Relief}: Relief Ratio, numeric
#' \item \strong{GrainSize}: Coarse Medium or Fine, 
#' \item \strong{Sample}: Field Sampling, A, B or C 
#' \item \strong{Points}: The number of points measured for each sample
#' \item \strong{Qm}: Multicrystalline Quarts Amount, percentage
#' \item \strong{Qp}: Polycrystalline Quarts Amount, percentage
#' \item \strong{Rf}: Rock Fragments Amount, percentage
#' \item \strong{M}: Mica Amount, percentage
#' }
#' @name data_sets_Fragments
#' @rdname data_sets_Fragments
#' @aliases Fragments fragments
#' @author Jeremy Hummon Grantham and Michael Anthony Velbel
#' @examples 
#' data(Fragments)
#' ggtern(Fragments,aes(Qm+Qp,Rf,M,colour=Sample)) +
#'   geom_density_tern(h=2,aes(fill=..level..),expand=0.75,alpha=0.5) + 
#'   geom_point(aes(shape=Position,size=Relief)) + 
#'   theme_bw() + theme_showarrows() + custom_percent('%') + 
#'   labs(title = "Grantham and Valbel Rock Fragment Data",
#'        x = "Q_m+Q_p", xarrow = "Quartz (Multi + Poly)",
#'        y = "R_f",     yarrow = "Rock Fragments",
#'        z = "M",       zarrow = "Mica") + 
#'   facet_wrap(~Sample,nrow=2)
NULL

