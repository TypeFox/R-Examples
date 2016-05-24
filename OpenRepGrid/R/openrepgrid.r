
############################### Package description ###########################

#' \pkg{OpenRepGrid}: an R package for the analysis of repertory grids. 
#'
#' The \pkg{OpenRepGrid} package provides tools for the analysis of repertory grid data.
#' The repertory grid is a method devised by George Alexander Kelly
#' in his seminal work "The Psychology of Personal Constructs" published in 1955.
#' The repertory grid has been used in and outside the context of Personal Construct 
#' Psychology (PCP) in a broad range of fields. For an introduction into the 
#' technique see e.g. Fransella, Bell and Bannister (2003).
#' 
#'
#' @author    Current members of the \pkg{OpenRepGrid} development team: Mark Heckmann. 
#'            Everyone who is interested in developing the package is invited to join.
#'
#'            The \pkg{OpenRepGrid} package development is hosted on github (\url{http://github.com/markheckmann/OpenRepGrid}).
#'            The github site provides information and allows to file bug reports or feature requests.
#'            Bug reports can also be emailed to the package maintainer or issued on 
#'            \url{http://www.openrepgrid.org} under section \emph{Suggestions/Issues}.
#'            The package maintainer is Mark Heckmann  <heckmann(at)uni-bremen.de>.
#'
#' @note      To get started with \pkg{OpenRepGrid} visit the project's home under \url{www.openrepgrid.org}. 
#'            On this site you will find tutorials, explanation about the theory, methods of analysis and 
#'            the according R code.
#'
#'            To see the preferable citation of the \pkg{OpenRepGrid} package, type 
#'            \code{citation("OpenRepGrid")} into the R console.
#'
#'            Disclaimer: Note that the package is distributed under the 
#'            \href{http://www.gnu.org/licenses/gpl-2.0.html}{GPL 2 license}.
#'            It is work in progress and is continuously being improved by hopefully 
#'            numerous contributors. It may contain bugs and errors.
#'            There is no warranty whatsoever for the program.
#'
#' @references  Fransella, F., Bell, R. C., & Bannister, D. (2003). \emph{A Manual for Repertory 
#'                  Grid Technique (2. Ed.)}. Chichester: John Wiley & Sons.
#'
#'              Kelly, G. A. (1955). \emph{The psychology of personal constructs. Vol. I, II.} 
#'                  New York: Norton, (2nd printing: 1991, Routledge, London, New York).
#' @keywords package repgrid
#' @name OpenRepGrid
#' @docType package
#' @import methods grid colorspace plyr stringr abind rgl GPArotation psych XML tcltk pvclust
NULL


#############################  Package overview  ##############################

#' \pkg{OpenRepGrid}: Annotated overview of package functions.   
#'
#' This documentation page contains an overview over the package functions
#' ordered by topics. The best place to start learning OpenRepGrid will
#' be the package website \url{http://www.openrepgrid.org} though.
#' 
#' @section Functions sorted by topic: 
#' 
#' \bold{Manipulating grids} \cr 
#' 
#' \tabular{ll}{
#'    \code{\link{left}}   \tab Move construct(s) to the left  \cr
#'    \code{\link{right}}  \tab Move construct(s) to the right \cr
#'    \code{\link{up}}     \tab Move construct(s) upwards \cr
#'    \code{\link{down}}   \tab Move construct(s) downwards \cr
#' }
#'
#' \bold{Loading and saving data} \cr 
#' 
#' \tabular{ll}{
#' \code{\link{importGridcor}}  \tab Import GRIDCOR data files \cr
#' \code{\link{importGridstat}}	\tab Import Gridstat data files \cr
#' \code{\link{importGridsuite}} \tab	Import Gridsuite data files \cr
#' \code{\link{importScivesco}}	\tab Import sci:vesco data files \cr
#' \code{\link{importTxt}}	\tab Import grid data from a text file \cr
#'  \tab \cr
#' \code{\link{saveAsTxt}} \tab Save grid in a text file (txt) \cr
#' }
#' 
#' \bold{Analyzing constructs} \cr 
#' 
#' Descriptive statistics of constructs
#' Construct correlations
#' distance
#' Root mean square of inter-construct correlations
#' Somers' D 
#' Principal component analysis (PCA) of construct correlation matrix 
#' Cluster analysis of constructs
#' 
#' \bold{Analyzing elements} \cr 
#' 
#' \bold{Visual representation} \cr 
#' 
#' \tabular{ll}{
#' \emph{Bertin plots} \tab \cr
#'   \tab \cr
#'   \code{\link{bertin}}             \tab  Make Bertin display of grid data \cr
#'   \code{\link{bertinCluster}}      \tab	Bertin display with corresponding cluster anaylsis \cr
#'   \tab \cr
#'   \tab \cr
#' \emph{Biplots} \tab \cr
#'   \tab \cr
#' \code{\link{biplot2d}}             \tab Draw a two-dimensional biplot \cr
#' \code{\link{biplotEsa2d}}          \tab Plot an eigenstructure analysis (ESA) biplot in 2D \cr
#' \code{\link{biplotSlater2d}}       \tab Draws Slater's INGRID biplot in 2D \cr
#'    \tab \cr
#' \code{\link{biplotPseudo3d}}       \tab See 'biplotPseudo3d' for its use. Draws a biplot of the grid in 2D with depth impression (pseudo 3D) \cr
#' \code{\link{biplotEsaPseudo3d}}    \tab Plot an eigenstructure analysis (ESA) in 2D grid with 3D impression (pseudo 3D) \cr
#' \code{\link{biplotSlaterPseudo3d}} \tab Draws Slater's biplot in 2D with depth impression (pseudo 3D) \cr
#'    \tab \cr
#' \code{\link{biplot3d}}	            \tab Draw grid in rgl (3D device) \cr
#' \code{\link{biplotEsa3d}}	        \tab Draw the eigenstructure analysis (ESA) biplot in rgl (3D device) \cr
#' \code{\link{biplotSlater3d}}	      \tab Draw the Slater's INGRID biplot in rgl (3D device) \cr
#'    \tab \cr
#' \code{\link{biplotSimple}}         \tab A graphically unsophisticated version of a biplot \cr
#' } 
#' 
#' \bold{Index measures} \cr
#'  
#' \tabular{ll}{
#' \code{\link{indexConflict1}}	  \tab Conflict measure for grids (Slade & Sheehan, 1979) based on correlations \cr
#' \code{\link{indexConflict2}}	  \tab Conflict measure for grids (Bassler et al., 1992) based on correlations \cr
#' \code{\link{indexConflict3}}	  \tab Conflict or inconsistenciy measure for grids (Bell, 2004) based on distances \cr
#' \code{\link{indexDilemma}}	    \tab Detect implicative dilemmas (conflicts) \cr
#'    \tab \cr
#' \code{\link{indexIntensity}}	  \tab Intensity index \cr
#' \code{\link{indexPvaff}}	      \tab Percentage of Variance Accounted for by the First Factor (PVAFF) \cr
#'    \tab \cr
#' \code{\link{indexBias}}        \tab Calculate 'bias' of grid as defined by Slater (1977) \cr
#' \code{\link{indexVariability}}	\tab Calculate 'variability' of a grid as defined by Slater (1977) \cr
#' }
#' 
#' \bold{Special features} \cr
#' 
#' \tabular{ll}{
#' \code{\link{alignByIdeal}}     \tab  Align constructs using the ideal element to gain pole preferences \cr
#' \code{\link{alignByLoadings}}  \tab	Align constructs by loadings on first pricipal component \cr
#' \code{\link{reorder2d}}        \tab Order grid by angles between construct and/or elements in 2D \cr
#' }
#' 
#' @section Settings:
#' 
#' \pkg{OpenRepGrid} uses several default settings e.g. to determine 
#' how many construct characters to display by default when displaying a grid.
#' The function \code{settings} can be used to show and change these settings.
#' Also it is possible to store the settings to a file and load the settings
#' file to restore the settings.
#' 
#' \tabular{ll}{
#' \code{\link{settings}}      \tab Show and modify global settings for OpenRepGrid \cr
#' \code{\link{settingsSave}}  \tab Save OpenRepGrid settings to file \cr
#' \code{\link{settingsLoad}}  \tab Load OpenRepGrid settings from file\cr
#' }
#' 
#' @section Grid datasets:
#' 
#' \pkg{OpenRepGrid} already contains some ready to use grid data sets. Most of 
#' the datasets are taken from the literature. To output the data simply type 
#' Type the name of the dataset to the console and press enter. \cr
#' 
#' \emph{Single grids} \cr
#' 
#' \tabular{ll}{
#' \code{\link{bell2010}}         \tab Grid data from a study by Haritos et al. (2004) 
#'                                     on role titles; used for demonstration of 
#'                                     construct alignment in Bell (2010, p. 46). \cr
#' \code{\link{bellmcgorry1992}}  \tab Grid from a psychotic patient used in Bell 
#'                                     (1997, p. 6). Data originated from a study 
#'                                     by Bell and McGorry (1992). \cr
#' \code{\link{boeker}}           \tab Grid from seventeen year old female schizophrenic 
#'                                     patient undergoing last stage of psychoanalytically 
#'                                     oriented psychotherapy (Boeker, 1996, p. 163). \cr
#' \code{\link{fbb2003}}          \tab Dataset used in \emph{A manual for Repertory Grid 
#'                                     Technique} (Fransella, Bell, & Bannister, 2003b, p. 60). \cr
#' \code{\link{feixas2004}}       \tab Grid from a 22 year old Spanish girl suffering 
#'                                     self-worth problems (Feixas & Saul, 2004, p. 77). \cr
#' \code{\link{mackay1992}}	      \tab Dataset \emph{Grid C} used in Mackay's paper on inter-element 
#'                                     correlation (1992, p. 65). \cr
#' \code{\link{leach2001a}}, \code{\link{leach2001b}} \tab	Pre- (a) and post-therapy (b) dataset from 
#'                                     sexual child abuse survivor (Leach, Freshwater, 
#'                                     Aldridge, & Sunderland, 2001, p. 227). \cr
#' \code{\link{raeithel}}         \tab Grid data to demonstrate the use of Bertin diagrams 
#'                                     (Raeithel, 1998, p. 223). The context of its 
#'                                     administration is unknown. \cr
#' \code{\link{slater1977a}}      \tab Drug addict grid dataset from (Slater, 1977, p. 32). \cr
#' \code{\link{slater1977b}}      \tab Grid dataset (ranked) from a seventeen year old 
#'                                     female psychiatric patient (Slater, 1977, p. 110) 
#'                                     showing depression, anxiety and self-mutilation. 
#'                                     The data was originally reported by Watson (1970).\cr
#' }
#' 
#' \emph{Multiple grids} \cr
#' 
#' NOT YET AVAILABLE \cr     
#' 
#' 
#' @section Functions for developers:
#' 
#' \pkg{OpenRepGrid}: internal functions overview for developers. \cr
#' 
#' Below you find a guide for developers: these functions are usually 
#' not needed by the casual user. The internal functions have a twofold goal
#' 1) to provide means for advanced numerical grid analysis and 2) 
#' to facilitate function development. The function for these purposes
#' are internal, i.e. they are not visible in the package documentation.
#' Nonetheless they do have a documentation that
#' can be accesses in the same way as for other functions.
#' More in the details section.
#' 
#' \bold{Functions for advanced grid analysis} \cr
#' 
#' The package provides functions to facilitate numerical research for grids. 
#' These comprise the generation of random data, permutation of grids etc. 
#' to facilitate Monte Carlo simulations, batch analysis of grids and other methods. 
#' With R as an underlying framework, the results of grid analysis easily lend 
#' themselves to further statistical processing and analysis within R. 
#' This is one of the central advantages for researchers compared to other 
#' standard grid software. The following table lists several functions for these purposes.
#'
#' \tabular{ll}{
#' \code{\link{randomGrid}}                       \tab  \cr
#' \code{\link{randomGrids}}                      \tab  \cr
#' \code{\link{permuteConstructs}}                \tab  \cr
#' \code{\link{permuteGrid}}                      \tab  \cr
#' \code{\link{quasiDistributionDistanceSlater}}  \tab  \cr
#' }
#'
#' \bold{Modules for function development} \cr
#' 
#' Beside the advanced analysis feature the developer's functions comprise 
#' low-level modules to create new functions for grid analysis. 
#' Though the internal structure of a repgrid object in R is simple 
#' (type e.g. \code{str(bell2010, 2)} to get an impression), it is convenient 
#' to not have to deal with access on this level. Several function like e.g. 
#' \code{getElementNames} are convenient wrappers that perform standard tasks 
#' needed when implementing new functions. The following table lists several 
#' functions for these purposes.
#'
#' \tabular{ll}{
#' \code{\link{getRatingLayer}}       \tab Retrieve grid scores from grid object. \cr
#' \code{\link{getNoOfConstructs}}    \tab Get the number of constructs in a grid object.    \cr
#' \code{\link{getNoOfElements}}      \tab Get the number of elements in a grid object.   \cr
#' \code{\link{dim}}                  \tab Get grid dimensions, i.e. contructs x elements.    \cr
#' \code{\link{getScale}}             \tab Get minimum and maximum scale value used in grid.  \cr
#' \code{\link{getScaleMidpoint}}     \tab Get midpoint of the grid rating scale.    \cr
#' \code{\link{getConstructNames}}    \tab Get construct names.                       \cr
#' \code{\link{getConstructNames2}}   \tab Get construct names (another newer version).      \cr
#' \code{\link{getElementNames}}      \tab Retrieve element names of repgrid object.  \cr
#' \code{\link{bindConstructs}}       \tab Concatenate the constructs of two grids.   \cr
#' \code{\link{doubleEntry}}          \tab  Join the constructs of a grid with the same reversed constructs.  \cr
#' }
#'
#' \bold{Other internal functions} \cr
#'
#' \tabular{ll}{
#' \code{\link{importTxtInternal}}  \tab  \cr
#' }
#'
#' @author    Current members of the \pkg{OpenRepGrid} development team: Mark Heckmann. 
#'            Everyone who is interested in developing the package is invited to join.
#'
#'            The \pkg{OpenRepGrid} package development is hosted on github (\url{http://github.com/markheckmann/OpenRepGrid}).
#'            The github site provides information and allows to file bug reports or feature requests.
#'            Bug reports can also be emailed to the package maintainer or issued on 
#'            \url{http://www.openrepgrid.org} under section \emph{Suggestions/Issues}.
#'            The package maintainer is Mark Heckmann <heckmann(at)uni-bremen.de>.
#'
#' @name OpenRepGrid-overview
#' @keywords package
#' @docType package
#'
NULL


