# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#
# 
#
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# .derive_flo:  generate the negative derivative for our florescence data;
#              this is done manually by measuring the delta between each
#              sequential data point
#
#   o Input:
#       -- flo.vec = numeric vector of florescence values
#   o Output:
#       -- numeric vector of derivative values
#          (length should be 331, but doesn't enforce)
#   o Usage:
#       -- .derive_flo (sample.flo.data)
#   o Example:
#       > sample.flo.data <- 1:131
#       > .derive_flo (sample.flo.data) [1:10]
#        [1] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
#
#   o Last modified:  06 May 2011
#   o Bugs:  None known
#   o Author:  C.A. Magaret
# --------------------------------------------------------------------------- #


#.derive_flo <- function (flo.vec) {
#    y.tmp.1 <- c (flo.vec[1], flo.vec)
#    y.tmp.2 <- c (flo.vec, flo.vec[length (flo.vec)])
#    y.tmp.1[2:length (y.tmp.1)] - y.tmp.2[2:length (y.tmp.2)]
#}
.derive_flo <- function (flo.yvec,flo.xvec) {
    diff.vec=vector(length=length(flo.yvec))
    for (i in 1:length(flo.yvec)-1)
    {
	diff.vec[i]= -(flo.yvec[i+1]-flo.yvec[i])/(flo.xvec[i+1]-flo.xvec[i])
    }
    diff.vec[length(flo.yvec)]=0
    return (diff.vec)
}

# --------------------------------------------------------------------------- #
# .scale_linear:  resclae the input vecotr to range between 0 and 1;
#
#   o Input:
#       -- inp.vec = numeric vector of values
#   o Output:
#       -- rescaled numeric vector of values
#          
#   o Usage:
#       -- .scale_linear (data)
#
#   o Example:
#
#   o Last modified:  11 May 2011
#   o Bugs:  None known
#   o Author:  D.A. Swan
# --------------------------------------------------------------------------- #


.scale_linear <- function (data.vec) {
    scale.ymin <- min (data.vec)
    scale.ymax <- max (data.vec)

    if (scale.ymin == scale.ymax)
	return (data.vec)

    scale.factor=scale.ymax-scale.ymin

    data2.vec=(data.vec-scale.ymin)/scale.factor
    return (data2.vec)
}

# --------------------------------------------------------------------------- #
# .load.flo / .load.abt:  read and parse the LightScanner output files
#
#   o Input:
#       -- file = character string of the path of the filename to be loaded
#   o Output:
#       -- data, as parsed from the indicated file, as matrix
#       -- (load.flo returns numeric matrix; load.abt returns character matrix)
#   o Usage:
#       -- load.flo (infile)
#       -- load.abt (infile)
#   o Example:
#       > abt.file.path <- "/data/file/belongs/here.ABT"
#       > flo.file.path <- "/data/file/belongs/here.FLO"
#       > data.abt <- load.abt (abt.file.path)
#       Read 480 items
#       > data.flo <- load.flo (flo.file.path)
#       Read 571968 items
#
#   o Last modified:  06 May 2011
#   o Bugs:  None known
#   o Author:  C.A. Magaret
# --------------------------------------------------------------------------- #


.load.abt <- function (file) {
    data <- matrix (scan (file=file, skip=7, what='character', nlines=96, sep="\t"), 
                    ncol=5, byrow=TRUE)
    abt.names <- c ('SampleName', 'WellNum', 'Unk1', 'Unk2', 'Unk3')
    names (data) <- abt.names
    return (data)
}


.load.flo <- function (file) {
    data <- matrix (as.numeric (scan (file=file, skip=1, what='numeric')), 
                    ncol=18, byrow=TRUE)
    flo.names <- c ('SampleNum', 'ProgramNum', 'SegmentNum', 'CycleNum', 'Time',
                    'Temperature', 'Flo1', 'Flo2', 'Flo3', 'Flo4', 'Flo5', 'Flo6', 
                    'Gain1', 'Gain2', 'Gain3', 'Gain4', 'Gain5', 'Gain6')
    names (data) <- flo.names
    return (data)
}

