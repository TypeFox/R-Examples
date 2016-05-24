##' emuR - Main Package of the EMU Speech Database Management System
##' 
##' The emuR package provides the next iteration of the EMU Speech 
##' Database Management System with database management, data 
##' extraction, data preparation and data visualization facilities. 
##' 
##' This package is part of the next iteration of the EMU Speech Database Management System (EMU_SDMS) 
##' which aims to be as close to an all-in-one solution for generating, manipulating, querying, 
##' analyzing and managing speech databases as possible. 
##' For an overview of the system please visit this URL: \url{http://ips-lmu.github.io/EMU.html}.
##'
##' It can be viewed as the main component of the EMU_SDMS as it acts as 
##' the central instance that is able to interact with every component of the system.
##' It takes care of database managing duties by being able to interact with a speech 
##' database that is stored in the emuDB format. Further, it has easy to understand and 
##' learn yet expressive and powerful querying mechanics, that allow the user to easily query 
##' the annotation structures of the database. Lastly it provides easy data extraction 
##' capabilities that extract data (e.g. formant values) which corresponds to the 
##' result of a query.
##' 
##' For an introduction to the emuR package please see the \code{emuR_intro} vignette 
##' by calling: \code{vignette('emuR_intro')} 
##' 
##' For information about the \code{emuDB} database format please see the \code{emuDB}
##' vignette by calling: \code{vignette('emuDB')}
##' 
##' For information about the query language used by the EMU_SDMS please see the \code{EQL}
##' vignette by calling: \code{vignette('EQL')}
##'    
##' Typical work-flow in emuR (emuDB required):
##' 
##' \enumerate{
##' \item Load database into current R session - \code{\link{load_emuDB}}
##' \item Database annotation / visual inspection - \code{\link{serve}} and connect the EMU-webApp to the local server
##' \item Query database - \code{\link{query}} (sometimes followed by \code{\link{requery_hier}} or \code{\link{requery_seq}})
##' \item Get trackdata (e.g. formant values) for the result of a query - \code{\link{get_trackdata}}
##' \item Data preparation
##' \item Visual data inspection
##' \item Further analysis and statistical processing
##' }
##' 
##' @name emuR-package
##' @aliases emuR emuR-package
##' 
##' @references Harrington, J. (2010). The Phonetic Analysis of Speech Corpora.
##' Blackwell.
##' 
##' @keywords package
##' @import methods
##' @docType package
##' @examples
##' \dontrun{
##' # create demo data including an emuDB called "ae" 
##' create_emuRdemoData(dir = tempdir())
##' 
##' # construct path to demo emuDB
##' path2ae = file.path(tempdir(), "emuR_demoData", "ae")
##' 
##' # load emuDB into current R session
##' ae = load_emuDB(path2ae)
##' 
##' # query loaded emuDB
##' lvowels = query(ae, "Phonetic = i: | u: | o:")
##' 
##' # extract labels from query result 
##' lvowels.labs = label(lvowels)
##' 
##' # list all ssffTrackDefinitions of emuDB
##' list_ssffTrackDefinitions(ae)
##' 
##' # get formant trackdata defined in ssffTrackDefinitions "fm" for query result
##' lvowels.fm = get_trackdata(ae, lvowels, "fm")
##' 
##' # extract track values at temporal midpoint of segments
##' lvowels.fmCut = dcut(lvowels.fm, .5, prop = TRUE)
##' 
##' # Plot the data as time signal and formant card
##' dplot(lvowels.fm[,1:2], lvowels.labs, normalise=TRUE, main = "Formants over vowel duration")
##' eplot(lvowels.fmCut[,1:2], lvowels.labs, dopoints=TRUE, 
##'       doellipse=FALSE, main = "F1/F2 of vowel midpoint", form=TRUE, 
##'       xlab = "F2 in Hz", ylab = "F1 in Hz")
##'       
##'       
##' # Plot of spectral data from 50% of aspiration duration
##' hs = query(ae,"Phonetic = H")
##' hs.labs = label(hs)
##' hs.dft = get_trackdata(ae, hs, "dft")
##' hs.dftCut = dcut(hs.dft, .5, prop=TRUE)
##' plot(hs.dftCut, hs.labs, main = "Spectral data of aspiration")
##' 
##' }
##' 
NULL



