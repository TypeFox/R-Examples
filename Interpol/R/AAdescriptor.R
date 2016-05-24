AAdescriptor <-
function(data,descriptor=151,normalize=0) {

	#################################################catch missing or bad input
	if(missing(data))
          stop("data argument missing")
        else
          if(!(is.character(data)))
            stop("invalid argument: data argument is required to be a vector of amino acids")
	if(!is.numeric(descriptor) | descriptor > 533 | descriptor < 1)
	  stop("invalid argument: descriptor is required to be 1 <= n <= 533")
	if(!is.numeric(normalize) | normalize > 2 | normalize < 0)
	  stop("invalid argument: normalization argument is required to be {1,2,3}")
	if(missing(descriptor))
		descriptor = 151
	if(missing(normalize))
		normalize = 0

	################################################Encoding
	
	#load("descriptors.RData")
	data(AAindex)
	indices <- AAindex[,1]
	AAindex = AAindex[,-1]

	descriptor = round(descriptor)
	normalize = round(normalize) 
	max_value = max(AAindex[descriptor,])

	################################################Determine normalization function to use
	if(normalize==0) { norm <- function(data,max_value) { return(data) } }
	if(normalize==1) { norm <- function(data,max_value) { return(data/max_value) } }
	if(normalize==2) { norm <- function(data,max_value) { return((data+max_value)/(max_value*2)) } }

	#data.encoded = list()
	data.encoded <- vector("list", length(data))
	for(n in 1:length(data)) {
		data.encoded.tmp = c()
		data.tmp = unlist(strsplit(data[n],""))
		for(m in 1:length(data.tmp)){ data.encoded.tmp = c(data.encoded.tmp, AAindex[descriptor, data.tmp[m]]) }
		data.encoded[[n]] <- norm(data.encoded.tmp, max_value)
	}
	names(data.encoded) <- names(data)

	return(data.encoded)
	
}

