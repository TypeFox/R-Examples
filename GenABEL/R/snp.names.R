"snp.names" <- 
function(data,begin,end,chromosome) {
	if (is(data,"gwaa.data") || is(data,"snp.data")) {
		out <- snp.names.gwaa.data(data,begin,end,chromosome)
		return(out)
	} else if (is(data,"scan.gwaa")) {
		if (missing(chromosome)) {
			if (missing(begin)) begin <- min(data$map)
			if (missing(end))   end <- max(data$map)
			if (is.character(begin)) begin <- data$map[data$snpname == begin]
			if (is.character(end)) end <- data$map[data$snpname == end]
			out <- data$snpname[data$map>=begin & data$map<=end]
			return(out)
		} else {
			if (missing(begin)) begin <- min(data$map[data$chromosome == chromosome])
			if (missing(end))   end <- max(data$map[data$chromosome == chromosome])
			if (is.character(begin)) begin <- data$map[data$snpname == begin]
			if (is.character(end)) end <- data$map[data$snpname == end]
			out <- data$snpname[data$map>=begin & data$map<=end & data$chromosome == chromosome]
			return(out)
		}
	} else if (is(data,"check.marker")) {
		if (missing(chromosome)) {
			if (missing(begin)) begin <- min(data$call$map)
			if (missing(end))   begin <- min(data$call$map)
			if (is.character(begin)) begin <- data$call$map[data$call$name == begin]
			if (is.character(end)) end <- data$call$map[data$call$name == end]
			out <- data$call$name[data$call$map>=begin & data$call$map<=end]
			return(out)
		} else {
			if (missing(begin)) begin <- min(data$call$map[data$call$chromosome == chromosome])
			if (missing(end))   begin <- min(data$call$map[data$call$chromosome == chromosome])
			if (is.character(begin)) begin <- data$call$map[data$call$name == begin]
			if (is.character(end)) end <- data$call$map[data$call$name == end]
			out <- data$call$name[data$call$map>=begin & data$call$map<=end & data$call$chromosome == chromosome]
			return(out)
		}
	} else {
		stop("data argument should be of type gwaa.data-class, snp.data-class, scan.gwaa-class or check.marker-class")
	}
}

"snp.names.gwaa.data" <- 
function(data,begin,end,chromosome) {
	if (is(data,"gwaa.data")) 
		data <- data@gtdata
	else if (!is(data,"snp.data")) stop("data should have class gwaa.data-class or snp.data-class")
	if (missing(begin) && missing(end)) {begin <- min(data@map); end <- max(data@map)}
	if (!missing(begin) && missing(end)) end <- max(data@map)
	if (missing(begin) && !missing(end)) begin <- min(data@map)
	if (!missing(chromosome)) {
		chromosome <- as.character(chromosome)
		if (length(chromosome)>1) stop("chromosome must be single number!")
		if (!any(data@chromosome == chromosome)) return(NULL)
		data <- data[,data@snpnames[data@chromosome == chromosome]]
	}
	if (is.character(begin)) {
		b1 <- data@map[begin]
		if (is.na(b1)) stop(paste("can not find marker",begin))
		begin <- b1
	}
	if (is.character(end)) {
		b1 <- data@map[end]
		if (is.na(b1)) stop(paste("can not find marker",end))
		end <- b1
	}
	if (exists("begin")) {
		if (length(begin)>1) stop("begin must be single number!")
		if (begin > max(data@map)) return(NULL)
	}
	if (exists("end")) {
		if (length(end)>1) stop("end must be single number!")
		if (end < min(data@map)) return(NULL)
	}
	if (begin > end) stop("Start-map parameter begin must be more than end-parameter")
	data <- data[,data@snpnames[data@map>=begin & data@map<=end]]
	out <- data@snpnames
	out
}

