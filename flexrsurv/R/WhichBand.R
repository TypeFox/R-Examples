# find the index of the band in which T is ( bands[whichband]< T <= bands[whichband+1]
# if length(bands) == Nbands,   0 <= whichband <= Nbands
# whichband == 0 if T <= bands[1]
# whichband == Nbands if T > bands[Nbands]

setGeneric("WhichBand",function(x, bands,...)standardGeneric("WhichBand"))
setMethod("WhichBand",
          signature("numeric", "GLMStepParam"),
          function(x, bands,...) bands@ncuts - findInterval(-x, -bands@cuts[bands@ncuts:1])
          )



setMethod("WhichBand",
          signature("numeric","numeric"),
          function(x, bands, ...)  length(bands) - findInterval(-x, -bands[length(bands):1])
          )

# find the index of the band in which T is ( bands[whichband]<= T < bands[whichband+1]
# if length(bands) == Nbands,   0 <= whichband <= Nbands
# whichband == 0 if T < bands[1]
# whichband == Nbands if T >= bands[Nbands]

setGeneric("WhichBandInf",function(x, bands,...)standardGeneric("WhichBandInf"))
setMethod("WhichBandInf",
          signature("numeric", "GLMStepParam"),
          function(x, bands,...) findInterval(x, bands@cuts)
          )



setMethod("WhichBandInf",
          signature("numeric","numeric"),
          function(x, bands, ...)  findInterval(x, bands)
          )

