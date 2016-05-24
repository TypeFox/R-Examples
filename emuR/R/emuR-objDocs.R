##' Segment list
##' 
##' A segment list is the result type of legacy Emu query.
##' 
##' 
##' @aliases segmentlist emusegs
##' @format multi-columned matrix one row per segment 
##' \itemize{ 
##'   \item columnlabel 
##'   \item columnsegment onset time 
##'   \item columnsegment offset time 
##'   \item columnutterance name 
##' }
##' @seealso \code{\link{query}}, \code{\link{demo.vowels}}
##' @keywords classes
##' @name segmentlist
##' @examples
##' 
##'    data(demo.vowels)
##'    
##'    #demo.vowels is a segment list
##'    demo.vowels
##' 
NULL

##' emuR segment list
##' @description
##' An emuR segment list is a list of segment descriptors. Each segment descriptor describes a sequence of annotation elements. The list is usually a result of an emuDB query using function \code{\link{query}}.
##' 
##' @details
##' Each row shows the annotation label sequence, the start and end position in time, session and bundle names, level name and type.
##' Additionally the row contains the UUID of the emuDB, the ID's of start and end elements and the corresponding start and end position as sample count and the sample rate.
##' These columns are not printed by default. The print method of emuRsegs hides them. To print all columns of a segment list object use the print method of \code{\link{data.frame}}.
##' For example to print all columns of an emuRsegs segmentlist \code{sl} type:
##' \code{print.data.frame(sl)}
##' Though the segment descriptors have references to the annotations, the label and sample/time position information is not updated if any of them change. The values of the segment list may get invalid if the the database is modified.
##' A segment may consist only of one single element, in this case start and end ID are equal.
##' An emuR segment list is the default result of \code{\link{query}} and can be used to get track data using \code{\link{get_trackdata}}.
##' The emuRsegs class inherits \link{emusegs} and hence \code{\link{data.frame}}
##' 
##' @aliases segment list emuRsegs
##'
##' @format Attributed data.frame, one row per segment descriptor.
##' 
##' Data frame columns are:
##' \itemize{ 
##'   \item labels: sequenced labels of segment concatenated by '->'
##'   \item start: onset time in milliseconds
##'   \item end: offset time in milliseconds
##'   \item session: session name
##'   \item bundle: bundle name
##'   \item level: level name
##'   \item type: type of "segment" row: 'ITEM': symbolic item, 'EVENT': event item, 'SEGMENT': segment
##'
##' }
##' Additional hidden columns:
##' \itemize{
##'  \item utts: utterance name (for compatibility to \link{emusegs} class)
##'  \item db_uuid: UUID of emuDB
##'  \item startItemID: item ID of first element of sequence
##'  \item endItemID: item ID of last element of sequence
##'  \item sampleStart: start sample position
##'  \item sampleEnd: end sample position
##'  \item sampleRate: sample rate
##' }
##' 
##' Attributes:
##' \itemize{
##'   \item database: name of emuDB
##'   \item query: Query string
##'   \item type: type ('segment' or 'event') (for compatibility to \link{emusegs} class)
##' }
##' 
##'
##' 
##' @seealso \code{\link{query}},\code{\link{get_trackdata}},\link{emusegs}
##' @keywords classes
##' @name emuRsegs
##' 
NULL





##' Start and end times for EMU segment lists and trackdata objects
##' 
##' Obtain start and end times for EMU segment lists and trackdata objects
##' 
##' The function returns the start and/or end times of either a segment list or
##' a trackdata object. The former refers to the boundary times of segments,
##' the latter the start and end times at which the tracks from segments occur.
##' start.emusegs and end.emusegs give exactly the same output as start and end
##' respectively.
##' 
##' @aliases start.emusegs end.emusegs start.trackdata end.trackdata
##' @param x a segment list or a trackdata object
##' @param ...  due to the generic only
##' @return A vector of times.
##' @author Jonathan Harrington
##' @seealso \code{\link{tracktimes}}
##' @keywords utilities
##' @name start.emusegs
##' @examples
##' 
##' # start time of a segment list
##' start(polhom)
##' # duration of a segment list
##' end(polhom) - start(polhom)
##' # duration from start time of segment list
##' # and start time of parallel EPG trackdata
##' start(polhom) - start(polhom.epg)
##' 
##' 
NULL





##' Track data object
##' 
##' A track data object is the result of get_trackdata().
##' 
##' 
##' @aliases trackdata Math.trackdata Math2.trackdata Ops.trackdata
##' Summary.trackdata
##' @format \describe{ \item{\$index}{a two columned matrix, each row keeps the
##' first and last index of the \$data rows that belong to one segment}
##' \item{\$ftime}{a two columned matrix, each row keeps the times marks of one
##' segment} \item{\$data}{a multi-columned matrix with the real track values
##' for each segment} }
##' @note The entire data track is retrieved for each segment in the segment
##' list. The amount of data returned will depend on the sample rate and number
##' of columns in the track requested.
##' @section Methods: The following generic methods are implemented for
##' trackdata obects.  \describe{ \item{list("Arith")}{\code{"+"}, \code{"-"},
##' \code{"*"}, \code{"^"}, \code{"\%\%"}, \code{"\%/\%"}, \code{"/"}}
##' \item{list("Compare")}{\code{"=="}, \code{">"}, \code{"<"}, \code{"!="},
##' \code{"<="}, \code{">="}} \item{list("Logic")}{\code{"&"}, \code{"|"}.  }
##' \item{list("Ops")}{\code{"Arith"}, \code{"Compare"}, \code{"Logic"}}
##' \item{list("Math")}{\code{"abs"}, \code{"sign"}, \code{"sqrt"},
##' \code{"ceiling"}, \code{"floor"}, \code{"trunc"}, X \code{"cummax"},
##' \code{"cummin"}, \code{"cumprod"}, \code{"cumsum"}, \code{"log"},
##' \code{"log10"}, \code{"log2"}, \code{"log1p"}, \code{"acos"},
##' \code{"acosh"}, \code{"asin"}, \code{"asinh"}, \code{"atan"},
##' \code{"atanh"}, \code{"exp"}, \code{"expm1"}, \code{"cos"}, \code{"cosh"},
##' \code{"sin"}, \code{"sinh"}, \code{"tan"}, \code{"tanh"}, \code{"gamma"},
##' \code{"lgamma"}, \code{"digamma"}, \code{"trigamma"} }
##' \item{list("Math2")}{\code{"round"}, \code{"signif"}}
##' \item{list("Summary")}{\code{"max"}, \code{"min"}, \code{"range"},
##' \code{"prod"}, \code{"sum"}, \code{"any"}, \code{"all"}} }
##' @seealso \code{\link{get_trackdata}}, \code{\link{demo.vowels.fm}}
##' \code{\link{demo.all.rms}}
##' @keywords classes
##' @name trackdata
##' @examples
##' 
##'    data(demo.vowels.fm)
##'    data(demo.vowels)
##'    
##'    #Formant track data for the first segment of the segment list demo.vowels
##'    demo.vowels.fm[1]
##'   
##' 
NULL

##' emuR track data object
##' 
##' A emuR track data object is the result of \code{\link{get_trackdata}} if the 
##' \code{resultType} parameter is set to \code{"emuRtrackdata"} or the result of 
##' an explicit call to \code{\link{create_emuRtrackdata}}. Compared to 
##' the \code{\link{trackdata}} object it is a sub-class of a \code{\link{data.table}}
##' / \code{\link{data.frame}} which is meant to ease integration with other
##' packages for further processing. It can be viewed as an amalgamation of
##' a \code{\link{emuRsegs}} and a \code{\link{trackdata}} object as it
##' contains the information stored in both objects.
##' 
##' 
##' @format The \code{\link{data.table}} / \code{\link{data.frame}} has the following columns:
##' 
##' \describe{
##'   \item{$sl_rowIdx}{column to indicate \code{\link{emuRsegs}} row index that 
##'                     the value belongs to}
##'   \item{$labels - $sampleRate}{duplicated information of \code{\link{emuRsegs}} row entries}
##'   \item{$times_rel}{relative time stamps of sample values in milliseconds}
##'   \item{$times_orig}{absolute time stamps of sample values in milliseconds}
##'   \item{$T1 - $TN}{actual data values (e.g. formant values / F0 values / DFT values / ...)}
##' }
##' 
##' Note that $labels - $sampleRate as well as $T1 - $TN (where the N in TN is to be read as the n-th T value) 
##' refer to multiple columns of the object.
##' 
##' @section Methods: The following methods are implemented for emuRtrackdata objects: 
##' 
##' \describe{ 
##'   \item{cut}{Function to extract a \code{\link{emuRtrackdata}} object from an emuRtrackdata at a single time point 
##'             or between two times}
##' }
##' @seealso \code{\link{get_trackdata}}, \code{\link{create_emuRtrackdata}}
##' @keywords classes
##' @name emuRtrackdata
##' @seealso trackdata
NULL
