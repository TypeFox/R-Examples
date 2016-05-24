.extractIndexes <-
function(nvarGroup){

	nvarGroupCumul <- cumsum(nvarGroup)
	f <- function(i){
		if(i==1){
			return(1:nvarGroupCumul[i])
		}else{
			return((nvarGroupCumul[i-1]+1):nvarGroupCumul[i])
		}
	}
	idx <- lapply(1:length(nvarGroupCumul), f)
	names(idx) <- names(nvarGroup)

	if(!all(as.numeric(summary(idx)[,"Length"]) == nvarGroup))
		stop("Extraction of indexes")
		idx
}
.extractIndexesSub <-
function(vect, idx){

	if(!is.character(vect[1]))
		stop("vect must be a vector of strings")
	if(!is.list(idx))
		stop("idx must be a list of indexes")

	f <- function(z){
		ii <- idx[[z]]
		if(is.null(ii))
			stop("wrong var names supplied")
		ii
	}
	lapply(vect, FUN=f)
}
.extractMember <-
function(obj, member){
	if(class(obj)!="list" | class(member)!="character")
		stop("Wrong classes")
	sapply(obj, FUN=function(l) l[[member]])
}


.reconstructSignal <- 
function(wavelet="d8", data){
	if(class(data)!="list")
		stop("the wavelet coefficients data must be a list")

	N <- length(unlist(data))
	series <- create.signalSeries()
	n.levels <- log2(N)
	filters <- wavDaubechies(wavelet = wavelet, normalized = FALSE)
	dict <- wavDictionary(wavelet = wavelet, dual = FALSE, decimate = FALSE, 
        n.sample = N, attr.x = NULL, n.levels = as.integer(n.levels), 
        boundary = "periodic", conv = TRUE, filters = filters, 
        fast = TRUE, is.complex = FALSE)

	waveletObject <- wavTransform(data = data, series = series, n.levels = as.integer(n.levels), 
        dictionary = dict, shifted = FALSE, xform = "dwt")
	reconstruct(waveletObject)
}

.vectorizeWavelets <- 
function(data){
	# l <- lapply(data[names(data)!="extra"], FUN=function(v) rev(v))
	l <- lapply(data, FUN=function(v) rev(v))
	rev(unlist(l))
}

