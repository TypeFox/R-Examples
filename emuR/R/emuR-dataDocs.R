
##' Three-columned matrix
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format Three-columned matrix
##' @name bridge
NULL





##' EPG-compressed trackdata from the segment list coutts
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name coutts.epg
NULL





##' Vector of word label from the segment list coutts
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name coutts.l
NULL





##' Segment list of words, read speech, female speaker of Australian English
##' from database epgcoutts
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name coutts
NULL





##' rms Data to coutts segment list
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @name coutts.rms
##' @format segmentlist
##' @examples
##' 
##' data(coutts.rms)
##' 
NULL





##' Trackdata of acoustic waveforms from the segment list coutts
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name coutts.sam
##' 
NULL





##' EPG-compressed trackdata from the segment list coutts2
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name coutts2.epg
NULL





##' Vector of word label from the segment list coutts2
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of word label
##' @name coutts2.l
NULL





##' Segment list, same as coutts but at a slower speech rate
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name coutts2
NULL





##' Trackdata of acoustic waveforms from the segment list coutts2
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name coutts2.sam
NULL





##' Emu segment list
##' 
##' Segment list of the demo database that is part of the Emu system.  It is
##' the result of a database query, that searched all segments at level
##' Phonetic.
##' 
##' A segment list is created via \code{\link{query}}.
##' 
##' @format First Column labels Second start time of the segment Third end time
##' of the segment Fourth utterance name of the utterance the segment was found
##' @seealso \code{\link{demo.vowels}} \code{\link{segmentlist}}
##' @docType data
##' @keywords datasets
##' @name demo.all
NULL





##' Emu track data for a rms track for segment list demo.all
##' 
##' A track list of the demo database that is part of the Emu system.  It is
##' the result of get rms data for the segment list demo.all (data(demo.all)).
##' 
##' A track list is created via the \code{\link{get_trackdata}} function.
##' 
##' @format A object with \$index, \$ftime and \$data
##' 
##' index: a two columned matrix with the range of the \$data rows that belong
##' to the segment ftime: a two columned matrix with the times marks of the
##' segment data: a vector with the rms data
##' @seealso \code{\link{demo.vowels.fm}} \code{\link{segmentlist}}
##' \code{\link{trackdata}}
##' @docType data
##' @keywords datasets
##' @name demo.all.rms
NULL





##' F0 track data for segment list demo.vowels
##' 
##' A track list of the demo database that is part of the Emu system.  It is
##' the result of get F0 data for the segment list demo.vowels (see
##' data(demo.vowels)).
##' 
##' A track list is created via the \code{\link{get_trackdata}} function.
##' 
##' @format An object with \$index, \$ftime and \$data
##' 
##' index: a two columned matrix with the range of the \$data rows that belong
##' to the segment ftime: a two columned matrix with the times marks of the
##' segment data: a one columned matrix with the F0 values
##' @seealso \code{\link{demo.all.rms}} \code{\link{segmentlist}}
##' \code{\link{trackdata}}
##' @docType data
##' @keywords datasets
##' @name demo.all.f0
NULL





##' Formant track data for segment list demo.vowels
##' 
##' A track list of the demo database that is part of the Emu system.  It is
##' the result of get fm data for the segment list demo.vowels (see
##' data(demo.vowels)).
##' 
##' A track list is created via the \code{\link{get_trackdata}} function.
##' 
##' @format index: a two columned matrix with the range of the \$data rows that
##' belong to the segment ftime: a two columned matrix with the times marks of
##' the segment data: a three columned matrix with the formant values of the
##' first three formants for each segment
##' @seealso \code{\link{demo.all.rms}} \code{\link{segmentlist}}
##' \code{\link{trackdata}}
##' @docType data
##' @keywords datasets
##' @name demo.all.fm
NULL





##' Emu segment List
##' 
##' Segment list of the demo database that is part of the Emu system.  It is
##' the result of a database query, that searched all vowel segments at level
##' Phonetic.
##' 
##' A segment list is created via \code{\link{query}}.
##' 
##' @format First Column labels Second start time of the segment Third end time
##' of the segment Fourth utterance name of the utterance the segment was found
##' @seealso \code{\link{demo.all}} \code{\link{segmentlist}}
##' @docType data
##' @keywords datasets
##' @name demo.vowels
NULL


#' F0 track data for segment list demo.vowels
#' 
#' A track list of the demo database that is part of the Emu system.  It is the
#' result of get F0 data for the segment list demo.vowels (see
#' data(demo.vowels)).
#' 
#' A track list is created via the \code{\link{get_trackdata}} function.
#' 
#' @format An object with \$index, \$ftime and \$data
#' 
#' index: a two columned matrix with the range of the \$data rows that belong
#' to the segment ftime: a two columned matrix with the times marks of the
#' segment data: a one columned matrix with the F0 values
#' @seealso \code{\link{demo.all.rms}} \code{\link{segmentlist}}
#' \code{\link{trackdata}}
#' @keywords datasets
#' @name demo.vowels.f0
NULL


#' Formant track data for segment list demo.vowels
#' 
#' A track list of the demo database that is part of the Emu system.  It is the
#' result of get fm data for the segment list demo.vowels (see
#' data(demo.vowels)).
#' 
#' A track list is created via the \code{\link{get_trackdata}} funciton.
#' 
#' @format index: a two columned matrix with the range of the \$data rows that
#' belong to the segment ftime: a two columned matrix with the times marks of
#' the segment data: a three columned matrix with the formant values of the
#' first three formants for each segment
#' @seealso \code{\link{demo.all.rms}} \code{\link{segmentlist}}
#' \code{\link{trackdata}}
#' @keywords datasets
#' @name demo.vowels.fm
NULL



##' Trackdata of formants from the segment list dip
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name dip.fdat
NULL





##' Vector of phoneme labels from the segment list dip
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of phoneme lables
##' @name dip.l
NULL





##' Segment list of dipththongs, two speakers one male, one female , Standard
##' North German, read speech from database kielread
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name dip
NULL





##' Vector of speaker labels from the segment list dip
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of speaker labels
##' @name dip.spkr
NULL







##' Spectral vector of a single E vowel produced by a male speaker of Standard
##' North German.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format spectral vector
##' @name e.dft
NULL





##' EPG-compressed trackdata from the segment list engassim
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name engassim.epg
NULL





##' Vector of phonetic labels from the segment list engassim: nK = nk,ng , sK =
##' sk,sg
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of phonetic labels
##' @name engassim.l
NULL





##' Segment list of a sequence of syllable final n or N preceding k or g ,
##' isolated words single speaker, Australian English female from database
##' epgassim.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name engassim
NULL



##' Vector of word labels from the segment list engassim.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of word labels
##' @name engassim.w
NULL





##' Spectral trackdata object from the segment list fric.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name fric.dft
NULL





##' Vector of labels from the segment list fric
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of labels
##' @name fric.l
NULL





##' Segment list of word-medial s or z one male speaker of Standard North
##' German, read speech from database kielread.
##' 
##' An EMU dataset
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name fric
NULL





##' Vector of word labels from the segment list fric.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of word labels
##' @name fric.w
NULL





##' Trackdata of formants from the segment list isol
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name isol.fdat
NULL





##' Vector of vowel phoneme labels from the segment list isol
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of vowel phoneme labels
##' @name isol.l
NULL





##' Segment list of vowels in a d d context isolated word speech, one male
##' speaker of Australian English from database isolated.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name isol
NULL





##' EPG-compressed trackdata from the segment list polhom
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name polhom.epg
NULL





##' Vector of phonetic labels from the segment list polhom
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of phonetic labels
##' @name polhom.l
NULL





##' Segment list of four Polish homorganic fricatives from database epgpolish.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name polhom
NULL







##' Data frame of various parameters and labels from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format dataframe
##' @name vowlax.df
NULL





##' Spectral matrix centred at the temporal midpoint of the vowels from the
##' segment list vowlax.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format spectral matrix
##' @name vowlax.dft.5
NULL





##' Matrix of formant data extracted at the temporal midpoint from the segment
##' list vowlax.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format matrix of formant data
##' @name vowlax.fdat.5
NULL





##' Trackdata of formants from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name vowlax.fdat
NULL





##' Vector of fundamental frequency extracted at the temporal midpoint from the
##' segment list vowlax.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of fundamental frequency
##' @name vowlax.fund.5
NULL





##' Trackdata of fundamental frequency from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name vowlax.fund
NULL





##' Vector of phoneme labels from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of phoenme labels
##' @name vowlax.l
NULL





##' Vector of labels preceding the vowels from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of phoenme labels
##' @name vowlax.left
NULL





##' Segment list of four lax vowels, read speech, one male and one female
##' speaker of Standard North German from database kielread.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format segmentlist
##' @name vowlax
NULL





##' Vector of labels following the vowels from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of phoenme labels
##' @name vowlax.right
NULL





##' Vector of RMS energy values at the temporal midpoint extracted at the
##' temporal midpoint from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of RMS energy values
##' @name vowlax.rms.5
NULL





##' Trackdata of RMS energy from the segment list vowlax
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format trackdata object
##' @name vowlax.rms
NULL





##' Vector of speaker labels from the segment list vowlax.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of speaker labels
##' @name vowlax.spkr
NULL





##' Vector of word labels from the segment list vowlax.
##' 
##' An EMU dataset
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of word labels
##' @name vowlax.word
NULL





##' Vector of word labels from segment list wordlax
##' 
##' For wordlax (see data(vowlax))
##' 
##' 
##' @docType data
##' @keywords datasets
##' @format vector of word labels
##' @name wordlax.l
NULL
