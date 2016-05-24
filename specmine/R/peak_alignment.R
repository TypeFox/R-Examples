
# group peaks with one of two methods
# method: "own" - own algorithm, "metaboanalyst" - metaboanalyst algorithm
# samp.classes: column of metadata dataframe - needed in case of "metaboanalyst"
# step used in "own" algorithm
# returns dataset using standard structure
"group_peaks" = function(sample.list, type, method = "own", metadata = NULL, samp.classes= 1, 
                         description = "", label.x = NULL, label.values = NULL, step = 0.03) {
  if (type == "nmr-peaks"){
	mzwid = 0.03
	bw = 10
  } else if (type == "lcms-peaks"){
	mzwid = 0.25
	bw = 30
  } else if (type == "gcms-peaks"){
	mzwid = 0.25
	bw = 5
  } else stop(paste("The ", type," type is not supported!", sep = ""))
	
  if (method == "own"){
		merged.peaks = merge_eq_peaks_samplelist(sample.list)
		samples.df = get_all_intensities(merged.peaks)
		data = group_peaks_own(samples.df, step)
	} 
  else if (method == "metaboanalyst"){
		data = group_peaks_metaboanalyst(sample.list, metadata[,samp.classes], names(sample.list), mzwid = mzwid, bw = bw)
	}
	create_dataset(as.matrix(data), metadata = metadata, type = type, 
                 description = description, label.x = label.x, label.values = label.values)
}


# grouping peaks: our algorithm
"group_peaks_own" = function(samples.df, step = 0.03)
{
  freqs = as.numeric(rownames(samples.df))
  st_intervals = create_intervals(freqs, step)
  initialized = F
  for (cent in st_intervals)
  {
    rows = samples.df[freqs >= cent & freqs < cent + step,]
    if (nrow(rows)==1)
    {
      if (initialized) {
        res = rbind(res, rows)
      }
      else { 
        res = rows
        initialized = T
      }
    }
    else if (nrow(rows) > 1)
    {
      newpeak = c()
      for(i in 1:ncol(samples.df)) {
        if (sum(!is.na(rows[,i])) > 0) newpeak[i] = sum(rows[,i], na.rm=T)
        else newpeak[i] = NA
      }
      if (initialized) {
        res = rbind(res, newpeak)
      }
      else {
        res = data.frame(t(newpeak))
        colnames(res) = colnames(samples.df)
        initialized = T
      }
      fs = c()
      for(k in 1:nrow(rows))
      {
        f = as.numeric(rownames(rows)[k])
        ocs = sum(!is.na(rows[k,]))
        fs = c(fs, rep(f, each = ocs))
      }
      m = round(median(fs),2)
      rownames(res)[nrow(res)] = as.character(m)
    }
  }
  names(res) = names(samples.df)
  res
}

# auxiliary function for our own grouping algorithm

"create_intervals" = function(freqs, step = 0.03)
{
  st = freqs[1]
  res = c(st)
  index = 1
  while(index < length(freqs))
  {
    while(index <= length(freqs) & freqs[index] < st + step) index = index + 1
    if (index <= length(freqs)) {
      st = freqs[index]
      res = c(res, st)
    }
  }
  res
}



# Group peak list based on position using xcms algorithm (align peaks wrt rt and mz)
# NMR peaks change ppm -> mz and add dummy rt
# 2-col MS need to add dummy rt
# 3-col MS can be used directly
# default mzwid MS 0.25 m/z, NMR 0.03 ppm
# bw 30 for LCMS, 5 for GCMS
group_peaks_metaboanalyst<-function(peaklist, samp.classes, samp.names,  mzwid = 0.03, bw = 10, minfrac = 0.5, minsamp = 1, max = 50) {
    samples <- samp.names;
    classlabel <- samp.classes;
    classnames <- levels(classlabel)

    classlabel <- as.vector(unclass(classlabel))
    classnum <- integer(max(classlabel))
    for (i in seq(along = classnum)){
        classnum[i] <- sum(classlabel == i)
    }

    peakmatrix = create_metaboanalyst_mat(peaklist)
    porder <- order(peakmatrix[,"ppm"])
    peakmat <- peakmatrix[porder,,drop=F]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])

    minpeakmat <- min(classnum)/2

    mass <- seq(peakmat[1,"ppm"], peakmat[nrow(peakmat),"ppm"] + mzwid, by = mzwid/2)
    masspos <- findEqualGreaterM(peakmat[,"ppm"], mass)

    groupmat <- matrix(nrow = 512, ncol = 7 + length(classnum))
    groupindex <- vector("list", 512)

    endidx <- 0
    num <- 0
    gcount <- integer(length(classnum))
    for (i in seq(length = length(mass)-2)) {
        startidx <- masspos[i]
        endidx <- masspos[i+2]-1
        if (endidx - startidx + 1 < minpeakmat)
            next
        speakmat <- peakmat[startidx:endidx,,drop=FALSE]
        den <- density(speakmat[,"rt"], bw, from = retrange[1]-3*bw, to = retrange[2]+3*bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2+length(classnum))
        snum <- 0
        while (deny[maxy <- which.max(deny)] > maxden/20 && snum < max) {
            grange <- descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(speakmat[,"rt"] >= den$x[grange[1]] & speakmat[,"rt"] <= den$x[grange[2]])
            gnum <- classlabel[unique(speakmat[gidx,"sample"])]
            for (j in seq(along = gcount))
                gcount[j] <- sum(gnum == j)
            if (! any(gcount >= classnum*minfrac & gcount >= minsamp))
                next
            snum <- snum + 1
            num <- num + 1
            ### Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            groupmat[num, 1] <- median(speakmat[gidx, "ppm"])
            groupmat[num, 2:3] <- range(speakmat[gidx, "ppm"])
            groupmat[num, 4] <- median(speakmat[gidx, "rt"])
            groupmat[num, 5:6] <- range(speakmat[gidx, "rt"])
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7+seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(porder[(startidx:endidx)[gidx]])
        }
    }
    colnames(groupmat) <- c("ppmmed", "ppmmin", "ppmmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", classnames)

    groupmat <- groupmat[seq(length = num),]
    groupindex <- groupindex[seq(length = num)]

    # Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])
    uindex <- rectUnique(groupmat[,c("ppmmin","ppmmax","rtmin","rtmax"),drop=FALSE],
                         uorder)
	
	res = list()
	res$groups <- groupmat[uindex,]
	res$groupidx<- groupindex[uindex]
	res$peaks <- peakmatrix
	set_groups_metaboanalyst(res, samp.names, samp.classes, ncol(peaklist[[1]]))
}

set_groups_metaboanalyst<-function(peaks.result, samp.names, samp.classes, num.col) {
    groupmat <- peaks.result$groups;
    groupindex <- peaks.result$groupidx;

    sampnum <- seq(length = length(samp.names))
    intcol <- match("int", colnames(peaks.result$peaks))
    sampcol <- match("sample", colnames(peaks.result$peaks))

    # row is peak, col is sample
    values <- matrix(nrow = length(groupindex), ncol = length(sampnum))

    for (i in seq(along = groupindex)) {
       # for each group, replace multiple peaks from the same sample by their sum
       for(m in sampnum){
            samp.inx<-which(peaks.result$peaks[groupindex[[i]], sampcol]==m)
            if(length(samp.inx)>0){
                 values[i, m] <- sum(peaks.result$peaks[groupindex[[i]][samp.inx], intcol]);
            }else{
                 values[i, m] <- NA;
            }
        }
    }
	values = data.frame(values)
    colnames(values) <- samp.names;
    if (num.col == 2){
		rownames(values) <- paste(round(groupmat[,paste("ppm", "med", sep="")],5));
	} else {
		rownames(values) <- paste(round(groupmat[,paste("ppm", "med", sep="")],5), "/", round(groupmat[,paste("rt", "med", sep="")],2), sep="");
	}
	values
}

# auxiliary functions for grouping using metaboanalyst algorithm 

# create matrix as used by MetaboAnalyst
create_metaboanalyst_mat <- function(sample.list){
  mat = matrix();
  allmat = NULL;
  if (ncol(sample.list[[1]]) == 2){
	for (i in 1:length(sample.list)){
		mat = cbind(sample.list[[i]][,1],1000, sample.list[[i]][,2],i)
		allmat = rbind(allmat,mat);
	}
  } else {
	for (i in 1:length(sample.list)){
		mat = cbind(sample.list[[i]][,1],sample.list[[i]][,2], sample.list[[i]][,3],i)
		allmat = rbind(allmat,mat);
	}
  }
  colnames(allmat) = c("ppm","rt","int", "sample")
  allmat  
}


descendMin = function (y, istart = which.max(y)) 
{
    if (!is.double(y)) 
        y <- as.double(y)
    unlist(.C("DescendMin", y, length(y), as.integer(istart - 
        1), ilower = integer(1), iupper = integer(1), DUP = FALSE, 
        PACKAGE = "xcms")[4:5]) + 1
}

findEqualGreaterM = function(x, val){
    idx = 1
    index = integer(length(val))
    for (i in 1:length(val)){
        while(idx <= length(x) && x[idx] < val[i]){
            idx = idx + 1
        }
        index[i] = idx
    }
    index
}
descendMin = function(y, istart = which.max(y)){
    if (!is.double(y)) 
        y <- as.double(y)
    istart = as.integer(istart)
    i = istart
    while (i > 1){
	if (y[i-1] >= y[i])
	    break
	i = i-1
    }
    ilower = i
    i = istart
    while (i < length(y)){
	if (y[i+1] >= y[i])
	    break
        i = i+1
    }
    iupper = i
    c(ilower= ilower, iupper=iupper)
}

rectUnique = function(m, order = seq(length = nrow(m)), xdiff = 0, ydiff = 0){
    nr = nrow(m)
    if (!is.double(m))
	m = as.double(m)
    x1 = 0
    x2 = nr
    y1 = nr * 2
    y2 = nr * 3
    xdiff = as.double(xdiff)
    ydiff = as.double(ydiff)
    keep = logical(nrow(m))
    i = 1

    while (i <= nr){
	io = order[i]
	keep[io] = 1
	j = 1
   	while (j<i){
	    jo = order[j]
	    if (keep[jo] && 
		!(m[x1+io] - m[x2+jo] > xdiff || m[x1+jo] - m[x2+io] > xdiff ||
		  m[y1+io] - m[y2+jo] > ydiff || m[y1+jo] - m[y2+io] > ydiff)) {
		keep[io] = 0
		break	
	    }
	    j = j+1
	}
    	i = i+1
    }
}
