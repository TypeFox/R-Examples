#####################
# VOWELS v. 1.2-1
# Vowel Manipulation, Normalization, and Plotting Package for R
# Tyler Kendall, tsk@uoregon.edu, 2009-2014
# cf. NORM: http://ncslaap.lib.ncsu.edu/tools/norm/ & http://lvc.uoregon.edu/norm/
#     Info: http://blogs.uoregon.edu/vowels/
#####################

#
# USAGE: source("vowels.R")
# --or-- R CMD install vowels_1.2-1.tar.gz (from command line)
#  Then library(vowels)
#  Then call functions as needed
#
#
# In all cases the input "vowels" is a data.frame of the following:
#  vowels<-data.frame(spkr, vtype, context_n, f1s, f2s, f3s, f1gs, f2gs, f3gs)
#
# The output of the normalization functions is the same, but with
#  attr(vowels, "norm.method") == the normalization method
#  attr(vowels, "no.f3s") == TRUE if the output data frame does NOT have
#   columns for F3 values, NA if there are F3 columns

#
# GENERAL FUNCTIONS
#

load.vowels <- function(file = NA, type='Hertz') {
	if (is.na(file)) file<-file.choose()
	vowels <- read.table(file, header=TRUE, sep="\t", quote="", comment.char="", blank.lines.skip=TRUE, fill=TRUE)
	attr(vowels, "unit.type")<-type
	vowels
}

label.columns <- function(vowels) {
	col3 <- "Context"
	vnames<-vector()
	if (!(is.null(attributes(vowels)$mean.values) | is.null(attributes(vowels)$mean.values))) col3 <- "N"
	vnames<-c("Speaker", "Vowel", col3, "F1", "F2", "F3", "F1.gl", "F2.gl", "F3.gl")
	if (!is.null(attributes(vowels)$no.f3s) | (dim(vowels)[2] < 9)){
		vnames<-c("Speaker", "Vowel", col3, "F1", "F2", "F1.gl", "F2.gl")
	}
	if (!is.null(attributes(vowels)$unit.type)) {
		if (attributes(vowels)$unit.type=="Bark") {
			vnames<-c("Speaker", "Vowel", col3, "Z1", "Z2", "Z3", "Z1.gl", "Z2.gl", "Z3.gl")
			if (!is.null(attributes(vowels)$no.f3s) | (dim(vowels)[2] < 9)) {
				vnames<-c("Speaker", "Vowel", col3, "Z1", "Z2", "Z1.gl", "Z2.gl")
			}
		} else if (attributes(vowels)$unit.type=="ERB") {
			vnames<-c("Speaker", "Vowel", col3, "E1", "E2", "E3", "E1.gl", "E2.gl", "E3.gl")
			if (!is.null(attributes(vowels)$no.f3s) | (dim(vowels)[2] < 9)) {
				vnames<-c("Speaker", "Vowel", col3, "E1", "E2", "E1.gl", "E2.gl")
			}
		}
	}
	names(vowels)<-vnames
	vowels
}

scalevowels <- function (normed.vowels) {
  
  f3.plus <- 0
  if (is.null(attributes(normed.vowels)$no.f3s)) {
  	f3.plus <- 1
  	min.f3 <- min(c(normed.vowels[,6], normed.vowels[,9]), na.rm=TRUE)
  	max.f3 <- max(c(normed.vowels[,6], normed.vowels[,9]), na.rm=TRUE)
  }
  
  min.f1 <- min(c(normed.vowels[,4], normed.vowels[,6+f3.plus]), na.rm=TRUE)
  max.f1 <- max(c(normed.vowels[,4], normed.vowels[,6+f3.plus]), na.rm=TRUE)
  min.f2 <- min(c(normed.vowels[,5], normed.vowels[,7+f3.plus]), na.rm=TRUE)
  max.f2 <- max(c(normed.vowels[,5], normed.vowels[,7+f3.plus]), na.rm=TRUE)
  
  normed.vowels[,4]<-round((500*(normed.vowels[,4] - min.f1)/(max.f1 - min.f1)) + 250, 3)
  normed.vowels[,6+f3.plus]<-round((500*(normed.vowels[,6+f3.plus] - min.f1)/(max.f1 - min.f1)) + 250, 3)
  normed.vowels[,5]<-round((1400*(normed.vowels[,5] - min.f2)/(max.f2 - min.f2)) + 850, 3)
  normed.vowels[,7+f3.plus]<-round((1400*(normed.vowels[,7+f3.plus] - min.f2)/(max.f2 - min.f2)) + 850, 3)
  if (f3.plus>0) {
  	normed.vowels[,6]<-round((1200*(normed.vowels[,6] - min.f1)/(max.f3 - min.f3)) + 2000, 3)
  	normed.vowels[,9]<-round((1200*(normed.vowels[,9] - min.f1)/(max.f3 - min.f3)) + 2000, 3)
  }

  attr(normed.vowels, "scaled.values")<-TRUE
  normed.vowels
}

compute.means <- function(vowels, separate=FALSE, speaker=NA) {

	if (!separate & !is.na(speaker)) {
		vals<-vowels[vowels[,1]==speaker,]
	} else {
		vals<-vowels
	}
	speakers<-unique(as.character(vals[,1]))
	if (!separate & (length(speakers) > 1)) {
		speakers<-c("AllSpkrs")
	}
	return.vals<-NA
	no.f3s<-FALSE
	if (!is.null(attributes(vowels)$no.f3s)) { no.f3s <- TRUE }
	for (s in speakers) {
		if (s == "AllSpkrs") {
			use.vals<-vals
		} else {
			use.vals<-vals[vals[,1]==s,]
		}
		vowel.types<-unique(use.vals[,2])
		spkrs<-vector(length=length(vowel.types))
		n.cols<-(length(use.vals)-3)
		if (n.cols > 6) n.cols <- 6
		formants<-matrix(nrow=length(vowel.types), ncol=n.cols)
		num.of.tokens<-vector(length=length(vowel.types))
		for (i in 1:length(vowel.types)) {
			spkrs[i]<-s
			if (length(speakers)==1) { spkrs[i]<-speakers[1] 
			} else if (!is.na(speaker)) { spkrs[i]<-speaker }
			
			for (j in 4:length(use.vals)) {
				if (is.numeric(use.vals[,j])) {
					formants[i,j-3]<-round(mean(use.vals[use.vals[,2]==vowel.types[i], j], na.rm=TRUE), 3)
				}
			}

			num.of.tokens[i]<-length(use.vals[,2][use.vals[,2]==vowel.types[i]])
		}
		if (is.data.frame(return.vals)) {
			new.vals<-data.frame(spkrs, vowel.types, num.of.tokens, formants)
			names(new.vals)<-c("Speaker", "Vowel", "N", names(vowels[,4:length(vowels)]))
			if (dim(formants)[2]==4) no.f3s<-TRUE
			return.vals<-rbind(return.vals, new.vals)
		} else {
			return.vals<-data.frame(spkrs, vowel.types, num.of.tokens, formants)
			names(return.vals)<-c("Speaker", "Vowel", "N", names(vowels[,4:length(vowels)]))
			if (dim(formants)[2]==4) no.f3s<-TRUE
		}
	}
	if (no.f3s) { attr(return.vals, "no.f3s")<-TRUE }
	attr(return.vals, "norm.method")<-attr(vowels, "norm.method")
	attr(return.vals, "norm.variant")<-attr(vowels,"norm.variant")
	attr(return.vals, "unit.type")<-attr(vowels,"unit.type")
	attr(return.vals, "mean.values")<-TRUE
	#return.vals<-label.columns(return.vals)
	return.vals
}

compute.sds <- function(vowels, separate=FALSE, speaker=NA) {

	if (!separate & !is.na(speaker)) {
		vals<-vowels[vowels[,1]==speaker,]
	} else {
		vals<-vowels
	}
	speakers<-unique(as.character(vals[,1]))
	if (!separate) {
		speakers<-c("AllSpkrs")
	}
	return.vals<-NA
	no.f3s<-FALSE
	if (!is.null(attributes(vowels)$no.f3s)) { no.f3s <- TRUE }
	for (s in speakers) {
		if (s == "AllSpkrs") {
			use.vals<-vals
		} else {
			use.vals<-vals[vals[,1]==s,]
		}
		vowel.types<-unique(use.vals[,2])
		spkrs<-vector(length=length(vowel.types))
		formants<-matrix(nrow=length(vowel.types), ncol=(length(use.vals)-3))
		num.of.tokens<-vector(length=length(vowel.types))
		for (i in 1:length(vowel.types)) {
			spkrs[i]<-s
			if (length(speakers)==1) { spkrs[i]<-speakers[1] 
			} else if (!is.na(speaker)) { spkrs[i]<-speaker }
			
			for (j in 4:length(use.vals)) {
				if (is.numeric(use.vals[,j])) {
					formants[i,j-3]<-round(sd(use.vals[use.vals[,2]==vowel.types[i], j], na.rm=TRUE), 3)
				}
			}

			num.of.tokens[i]<-length(use.vals[,2][use.vals[,2]==vowel.types[i]])
		}
		if (is.data.frame(return.vals)) {
			new.vals<-data.frame(spkrs, vowel.types, num.of.tokens, formants)
			names(new.vals)<-c("Speaker", "Vowel", "N", names(vowels[,4:length(vowels)]))
			if (dim(formants)[2]==4) no.f3s<-TRUE
			return.vals<-rbind(return.vals, new.vals)
		} else {
			return.vals<-data.frame(spkrs, vowel.types, num.of.tokens, formants)
			names(return.vals)<-c("Speaker", "Vowel", "N", names(vowels[,4:length(vowels)]))
			if (dim(formants)[2]==4) no.f3s<-TRUE
		}
	}
	if (no.f3s) { attr(return.vals, "no.f3s")<-TRUE }
	attr(return.vals, "norm.method")<-attr(vowels, "norm.method")
	attr(return.vals, "norm.variant")<-attr(vowels,"norm.variant")
	attr(return.vals, "unit.type")<-attr(vowels,"unit.type")
	attr(return.vals, "standard.devs")<-TRUE
	return.vals
}

compute.medians <- function(vowels, separate=FALSE, speaker=NA) {

	if (!separate & !is.na(speaker)) {
		vals<-vowels[vowels[,1]==speaker,]
	} else {
		vals<-vowels
	}
	speakers<-unique(as.character(vals[,1]))
	if (!separate & (length(speakers) > 1)) {
		speakers<-c("AllSpkrs")
	}
	return.vals<-NA
	no.f3s<-FALSE
	if (!is.null(attributes(vowels)$no.f3s)) { no.f3s <- TRUE }
	for (s in speakers) {
		if (s == "AllSpkrs") {
			use.vals<-vals
		} else {
			use.vals<-vals[vals[,1]==s,]
		}
		vowel.types<-unique(use.vals[,2])
		spkrs<-vector(length=length(vowel.types))
		n.cols<-(length(use.vals)-3)
		if (n.cols > 6) n.cols <- 6
		formants<-matrix(nrow=length(vowel.types), ncol=n.cols)
		num.of.tokens<-vector(length=length(vowel.types))
		for (i in 1:length(vowel.types)) {
			spkrs[i]<-s
			if (length(speakers)==1) { spkrs[i]<-speakers[1] 
			} else if (!is.na(speaker)) { spkrs[i]<-speaker }
			
			for (j in 4:length(use.vals)) {
				if (is.numeric(use.vals[,j])) {
					formants[i,j-3]<-round(median(use.vals[use.vals[,2]==vowel.types[i], j], na.rm=TRUE), 3)
				}
			}

			num.of.tokens[i]<-length(use.vals[,2][use.vals[,2]==vowel.types[i]])
		}
		if (is.data.frame(return.vals)) {
			new.vals<-data.frame(spkrs, vowel.types, num.of.tokens, formants)
			names(new.vals)<-c("Speaker", "Vowel", "N", names(vowels[,4:length(vowels)]))
			if (dim(formants)[2]==4) no.f3s<-TRUE
			return.vals<-rbind(return.vals, new.vals)
		} else {
			return.vals<-data.frame(spkrs, vowel.types, num.of.tokens, formants)
			names(return.vals)<-c("Speaker", "Vowel", "N", names(vowels[,4:length(vowels)]))
			if (dim(formants)[2]==4) no.f3s<-TRUE
		}
	}
	if (no.f3s) { attr(return.vals, "no.f3s")<-TRUE }
	attr(return.vals, "norm.method")<-attr(vowels, "norm.method")
	attr(return.vals, "norm.variant")<-attr(vowels,"norm.variant")
	attr(return.vals, "unit.type")<-attr(vowels,"unit.type")
	attr(return.vals, "median.values")<-TRUE
	#return.vals<-label.columns(return.vals)
	return.vals
}

#
# CONVERSION FUNCTIONS
# Based on Traunmuller (1997)'s formula
# These conversion functions all take vowel dataframes in Hz

convert.bark <- function(vowels) {
	z1s<-round((26.81/(1+1960/vowels[,4]))-0.53,3)
	z2s<-round((26.81/(1+1960/vowels[,5]))-0.53,3)
	z3s<-round((26.81/(1+1960/vowels[,6]))-0.53,3)

	z1.gls<-round((26.81/(1+1960/vowels[,7]))-0.53,3)
	z2.gls<-round((26.81/(1+1960/vowels[,8]))-0.53,3)
	z3.gls<-round((26.81/(1+1960/vowels[,9]))-0.53,3)

	return.vals<-data.frame(vowels[,1], vowels[,2], vowels[,3], z1s, z2s, z3s, z1.gls, z2.gls, z3.gls)
	#names(return.vals)<-c("Speaker", "Vowel", "Context", "Z1", "Z2", "Z3", "Z1 gl", "Z2 gl", "Z3 gl")
	attr(return.vals, "unit.type")<-"Bark"
	attr(return.vals, "mean.values")<-attr(vowels,"mean.values")
	attr(return.vals, "median.values")<-attr(vowels,"median.values")
	return.vals<-label.columns(return.vals)
	return.vals
}

convert.erb <- function(vowels) {
	e1s<-round(11.17*log((vowels[,4]+312)/(vowels[,4]+14675))+43.0,3)
	e2s<-round(11.17*log((vowels[,5]+312)/(vowels[,5]+14675))+43.0,3)
	e3s<-round(11.17*log((vowels[,6]+312)/(vowels[,6]+14675))+43.0,3)

	e1.gls<-round(11.17*log((vowels[,7]+312)/(vowels[,7]+14675))+43.0,3)
	e2.gls<-round(11.17*log((vowels[,8]+312)/(vowels[,8]+14675))+43.0,3)
	e3.gls<-round(11.17*log((vowels[,9]+312)/(vowels[,9]+14675))+43.0,3)

	return.vals<-data.frame(vowels[,1], vowels[,2], vowels[,3], e1s, e2s, e3s, e1.gls, e2.gls, e3.gls)
	#names(return.vals)<-c("Speaker", "Vowel", "Context", "E1", "E2", "E3", "E1 gl", "E2 gl", "E3 gl")
	attr(return.vals, "unit.type")<-"ERB"
	attr(return.vals, "mean.values")<-attr(vowels,"mean.values")
	attr(return.vals, "median.values")<-attr(vowels,"median.values")
	return.vals<-label.columns(return.vals)
	return.vals
}

#
# NORMALIZATION FUNCTIONS
# 

# Thomas' variation of Syrdal & Gopal
norm.bark <- function(vowels) {
	# convert formant mean values to Bark
 	bvowels<-convert.bark(vowels)
 	if (!is.null(attributes(vowels)$unit.type)) {
	  	if (attributes(vowels)$unit.type=="Bark") {
  			bvowels<-vowels
  		}
  	}
	z1s<-bvowels[,4]
	z2s<-bvowels[,5]
	z3s<-bvowels[,6]

	z1.gls<-bvowels[,7]
	z2.gls<-bvowels[,8]
	z3.gls<-bvowels[,9]
	
	# determine z3-z1, z3-z2, and z2-z1
	z3.z1<-round(z3s - z1s, 3)
	z3.z2<-round(z3s - z2s, 3)
	z2.z1<-round(z2s - z1s, 3)
	
	z3.gl.z1<-round(z3.gls - z1.gls, 3)
	z3.gl.z2<-round(z3.gls - z2.gls, 3)
	z2.gl.z1<-round(z2.gls - z1.gls, 3)
	return.vals<-data.frame(vowels[,1], vowels[,2], vowels[,3], z3.z1, z3.z2, z2.z1, z3.gl.z1, z3.gl.z2, z2.gl.z1)
	names(return.vals)<-c("Speaker", "Vowel", "Context", "Z3-Z1", "Z3-Z2", "Z2-Z1", "Z3-Z1 gl", "Z3-Z2 gl", "Z2-Z1 gl")
	attr(return.vals, "norm.method")<-"Bark Difference"
	attr(return.vals, "unit.type")<-"Bark"
	return.vals
}

norm.labov <- function(vowels, G.value = NA, use.f3 = FALSE, geomean = TRUE) {
	
	if (is.na(G.value)) G.value <- 6.896874
	
	return.vals<-NA
	for (speaker in as.character(unique(vowels[,1]))) {
		vals<-vowels[vowels[,1]==speaker,]
	
		# do Labov's Nearey-like algorithm
		ln.f1s <- log(vals[,4])
		ln.f2s <- log(vals[,5])
		ln.f1.gls <- log(vals[,7])
		ln.f2.gls <- log(vals[,8])
		# set all NA values to be the same as onset, basically act  
		# like all monophthongs are diphs with glide same as onset 
		ln.f1.gls[is.na(ln.f1.gls)]<-ln.f1s[is.na(ln.f1.gls)]
		ln.f2.gls[is.na(ln.f2.gls)]<-ln.f2s[is.na(ln.f2.gls)]
		
		if (use.f3) {
			ln.f3s <- log(vals[,6])
			ln.f3.gls <- log(vals[,9])
			ln.f3.gls[is.na(ln.f3.gls)]<-ln.f3s[is.na(ln.f3.gls)]
		
			grand.mean <- mean(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls, ln.f3s, ln.f3.gls), na.rm=TRUE)
			if (geomean) grand.mean <- exp(mean(log(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls, ln.f3s, ln.f3.gls)), na.rm=TRUE))
		} else {
			grand.mean <- mean(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls), na.rm=TRUE)
			if (geomean) grand.mean <- exp(mean(log(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls)), na.rm=TRUE))
		}
	
	    factor <- exp(G.value - grand.mean)
		normed.f1s <- round(factor * vals[,4], 3)
		normed.f2s <- round(factor * vals[,5], 3)
		normed.f1.gls <- round(factor * vals[,7], 3)
		normed.f2.gls <- round(factor * vals[,8], 3)
		
		normed.f3s<-NA
		normed.f3.gls<-NA
		if (use.f3) {
			normed.f3s <- round(factor * vals[,6], 3)
			normed.f3.gls <- round(factor * vals[,9], 3)
		
			if (is.data.frame(return.vals)) {
				new.vals<-data.frame(vals[,1], vals[,2], vals[,3], normed.f1s, normed.f2s, normed.f3s, normed.f1.gls, normed.f2.gls, normed.f3.gls)
				names(new.vals)<-c("Speaker", "Vowel", "Context", "F1'", "F2'", "F3'", "F1' gl", "F2' gl", "F3' gl")
				return.vals<-rbind(return.vals, new.vals)
			} else {
				return.vals<-data.frame(vals[,1], vals[,2], vals[,3], normed.f1s, normed.f2s, normed.f3s, normed.f1.gls, normed.f2.gls, normed.f3.gls)
				names(return.vals)<-c("Speaker", "Vowel", "Context", "F1'", "F2'", "F3'", "F1' gl", "F2' gl", "F3' gl")
			}
		} else {
			if (is.data.frame(return.vals)) {
				new.vals<-data.frame(vals[,1], vals[,2], vals[,3], normed.f1s, normed.f2s, normed.f1.gls, normed.f2.gls)
				names(new.vals)<-c("Speaker", "Vowel", "Context", "F1'", "F2'", "F1' gl", "F2' gl")
				return.vals<-rbind(return.vals, new.vals)
			} else {
				return.vals<-data.frame(vals[,1], vals[,2], vals[,3], normed.f1s, normed.f2s, normed.f1.gls, normed.f2.gls)
				names(return.vals)<-c("Speaker", "Vowel", "Context", "F1'", "F2'", "F1' gl", "F2' gl")
			}
			attr(return.vals, "no.f3s")<-TRUE
		}
	}
	attr(return.vals, "norm.method")<-"Labov"
	if (G.value==6.896874) {
		attr(return.vals, "norm.variant")<-"Labov, w/ Telsur G"
	}
	attr(return.vals, "G.value")<-G.value
	return.vals
}

norm.lobanov <- function(vowels, f1.all.mean=NA, f2.all.mean=NA) {
	
	return.vals<-NA
	for (speaker in as.character(unique(vowels[,1]))) {
		vals<-vowels[vowels[,1]==speaker,]
	
		f1s<-vals[,4]
		f2s<-vals[,5]
		f1.gls<-vals[,7]
		f2.gls<-vals[,8]
		f1.gls[is.na(f1.gls)]<-f1s[is.na(f1.gls)]
		f2.gls[is.na(f2.gls)]<-f2s[is.na(f2.gls)]
		mean.f1 <- mean(c(f1s, f1.gls))
    	if (!is.na(f1.all.mean)) mean.f1<-f1.all.mean
		mean.f2 <- mean(c(f2s, f2.gls))
    	if (!is.na(f2.all.mean)) mean.f2<-f2.all.mean
		stdev.f1 <- sd(c(f1s, f1.gls))
		stdev.f2 <- sd(c(f2s, f2.gls))
		norm.f1s <- round((f1s-mean.f1)/stdev.f1, 3)
		norm.f2s <- round((f2s-mean.f2)/stdev.f2, 3)
		norm.f1.gls <- round((f1.gls-mean.f1)/stdev.f1, 3)
		norm.f2.gls <- round((f2.gls-mean.f2)/stdev.f2, 3)
		norm.f1.gls[is.na(vals[,7])]<-NA
		norm.f2.gls[is.na(vals[,8])]<-NA	
		
		if (is.data.frame(return.vals)) {
			new.vals<-data.frame(vals[,1], vals[,2], vals[,3], norm.f1s, norm.f2s, norm.f1.gls, norm.f2.gls)
			names(new.vals)<-c("Speaker", "Vowel", "Context", "F*1", "F*2", "F*1 gl", "F*2 gl")
			return.vals<-rbind(return.vals, new.vals)
		} else {
			return.vals<-data.frame(vals[,1], vals[,2], vals[,3], norm.f1s, norm.f2s, norm.f1.gls, norm.f2.gls)
			names(return.vals)<-c("Speaker", "Vowel", "Context", "F*1", "F*2", "F*1 gl", "F*2 gl")
		}
	}
	attr(return.vals, "norm.method")<-"Lobanov"
	attr(return.vals, "no.f3s")<-TRUE
	return.vals
}

norm.nearey <- function(vowels, formant.int = FALSE, use.f3 = FALSE, all.mean = NA, all.mean.f2 = NA, all.mean.f3 = NA) {
	
	return.vals<-NA
	for (speaker in as.character(unique(vowels[,1]))) {
		vals<-vowels[vowels[,1]==speaker,]

		# do Nearey algorithm on the value
		ln.f1s <- log(vals[,4])
		ln.f2s <- log(vals[,5])
		ln.f1.gls <- log(vals[,7])
		ln.f2.gls <- log(vals[,8])
		# set all NA values to be the same as onset, basically act like 
		# all monophthongs are diphs with glide same as onset 
		ln.f1.gls[is.na(ln.f1.gls)]<-ln.f1s[is.na(ln.f1.gls)]
		ln.f2.gls[is.na(ln.f2.gls)]<-ln.f2s[is.na(ln.f2.gls)]
		grand.f1 <- mean(c(ln.f1s, ln.f1.gls), na.rm=TRUE)
		grand.f2 <- mean(c(ln.f2s, ln.f2.gls), na.rm=TRUE)
		
		ln.f3s<-vector(length=length(ln.f1s))
		ln.f3.gls<-vector(length=length(ln.f1.gls))
		grand.f3 <- NA
		if (use.f3) {
			ln.f3s <- log(vals[,6])
			ln.f3.gls <- log(vals[,9])
			ln.f3.gls[is.na(ln.f3.gls)]<-ln.f3s[is.na(ln.f3.gls)]
			grand.f3 <- mean(c(ln.f3s, ln.f3.gls), na.rm=TRUE)
		}
		
		if (!formant.int) {
		   	if (!use.f3) {
		   		grand.f1 <- mean(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls), na.rm=TRUE)
		   		grand.f2 <- mean(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls), na.rm=TRUE)
		   	} else {
				grand.f1 <- mean(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls, ln.f3s, ln.f3.gls), na.rm=TRUE)
		   		grand.f2 <- mean(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls, ln.f3s, ln.f3.gls), na.rm=TRUE)
		   		grand.f3 <- mean(c(ln.f1s, ln.f1.gls, ln.f2s, ln.f2.gls, ln.f3s, ln.f3.gls), na.rm=TRUE)
		   	}
		}
    	if (!is.na(all.mean)) {
		   grand.f1 <- all.mean
		   if (!is.na(all.mean.f2)) {
			   grand.f2 <- all.mean.f2	   
		   } else {
			   grand.f2 <- all.mean
		   }
		   if (!is.na(all.mean.f3)) {
			   grand.f3 <- all.mean.f3	   
		   } else {
			   grand.f3 <- all.mean
		   }
		}
		exp.ln.f1s <- round(exp(ln.f1s - grand.f1), 3)
		exp.ln.f2s <- round(exp(ln.f2s - grand.f2), 3)
		exp.ln.f1.gls <- round(exp(ln.f1.gls - grand.f1), 3)
		exp.ln.f2.gls <- round(exp(ln.f2.gls - grand.f2), 3)

		# reset NA glides to NA
		exp.ln.f1.gls[is.na(vals[,7])]<-NA
		exp.ln.f2.gls[is.na(vals[,8])]<-NA
	
		if (use.f3) {
			exp.ln.f3s <- round(exp(ln.f3s - grand.f3), 3)
			exp.ln.f3.gls <- round(exp(ln.f3.gls - grand.f3), 3)

			# reset NA glides to NA
			exp.ln.f3.gls[is.na(vals[,9])]<-NA
				
			if (is.data.frame(return.vals)) {
				new.vals<-data.frame(vals[,1], vals[,2], vals[,3], exp.ln.f1s, exp.ln.f2s, exp.ln.f3s, exp.ln.f1.gls, exp.ln.f2.gls, exp.ln.f3.gls)
				names(new.vals)<-c("Speaker", "Vowel", "Context", "F*1", "F*2", "F*3", "F*1 gl", "F*2 gl", "F*3 gl")
				return.vals<-rbind(return.vals, new.vals)
			} else {
				return.vals<-data.frame(vals[,1], vals[,2], vals[,3], exp.ln.f1s, exp.ln.f2s, exp.ln.f3s, exp.ln.f1.gls, exp.ln.f2.gls, exp.ln.f3.gls)
				names(return.vals)<-c("Speaker", "Vowel", "Context", "F*1", "F*2", "F*3", "F*1 gl", "F*2 gl", "F*3 gl")
			}
		} else {	
			if (is.data.frame(return.vals)) {
				new.vals<-data.frame(vals[,1], vals[,2], vals[,3], exp.ln.f1s, exp.ln.f2s, exp.ln.f1.gls, exp.ln.f2.gls)
				names(new.vals)<-c("Speaker", "Vowel", "Context", "F*1", "F*2", "F*1 gl", "F*2 gl")
				return.vals<-rbind(return.vals, new.vals)
			} else {
				return.vals<-data.frame(vals[,1], vals[,2], vals[,3,], exp.ln.f1s, exp.ln.f2s, exp.ln.f1.gls, exp.ln.f2.gls)
				names(return.vals)<-c("Speaker", "Vowel", "Context", "F*1", "F*2", "F*1 gl", "F*2 gl")
			}
			attr(return.vals, "no.f3s")<-TRUE
		}
	}
	attr(return.vals, "norm.method")<-"Nearey"
	if (formant.int) {
		attr(return.vals, "norm.variant")<-"Nearey1"
	} else {
		attr(return.vals, "norm.variant")<-"Nearey2"
	}
	return.vals
}

norm.wattfabricius <- function(vowels, norm.means=FALSE, mod.WF=FALSE) {
        	
	return.vals<-NA
	for (speaker in as.character(unique(vowels[,1]))) {
		vals<-vowels[vowels[,1]==speaker,]
		
		# do Watt & Fabricius algorithm

		# get idealized min&max F1 and F2, according to Watt & Fabricius
	
		mean.vals <- compute.means(vals)
		col.plus <- 0
		if (dim(mean.vals)[2]>7) { col.plus<-1 }
	
		# unlike W&F (2003), we determine the corners of the vowel triangle
		#  by choosing the vowels automatically based on position...
		i.f1 <- min(c(mean.vals[,4], mean.vals[,6+col.plus]), na.rm=TRUE)
		i.f2 <- max(c(mean.vals[,5], mean.vals[,7+col.plus]), na.rm=TRUE)
		ul.f1 <- i.f1
		ul.f2 <- i.f1
		a.f1 <- max(c(mean.vals[,4], mean.vals[,6+col.plus]), na.rm=TRUE)
		a.f2 <- mean(mean.vals[,5][mean.vals[,4]==a.f1 | mean.vals[,6+col.plus]==a.f1], na.rm=TRUE)
		
		if (mod.WF) { a.f2 <- (i.f2 + ul.f2)/2 }

		# compute S values
		S.f1 <- (i.f1 + a.f1 + ul.f1)/3
		S.f2 <- (i.f2 + a.f2 + ul.f2)/3
	
		# divide f1 and f2 values by S values
		if (norm.means) {
			norm.f1s <- round(mean.vals[,4]/S.f1, 3)
			norm.f2s <- round(mean.vals[,5]/S.f2, 3)
			norm.f1.gls <- round(mean.vals[,6+col.plus]/S.f1, 3)
			norm.f2.gls <- round(mean.vals[,7+col.plus]/S.f2, 3)
			
			if (is.data.frame(return.vals)) {
				new.vals <- data.frame(mean.vals[,1], mean.vals[,2], mean.vals[,3], norm.f1s, norm.f2s, norm.f1.gls, norm.f2.gls)
				names(new.vals) <- c("Speaker", "Vowel", "N", "F1/S(F1)", "F2/S(F2)", "F1g/S(F1g)", "F2g/S(F2g)")
				return.vals<-rbind(return.vals, new.vals)
			} else {
				return.vals <- data.frame(mean.vals[,1], mean.vals[,2], mean.vals[,3], norm.f1s, norm.f2s, norm.f1.gls, norm.f2.gls)
				names(return.vals) <- c("Speaker", "Vowel", "N", "F1/S(F1)", "F2/S(F2)", "F1g/S(F1g)", "F2g/S(F2g)")
			}
			attr(return.vals, "mean.values")<-TRUE
		} else {
			norm.f1s <- round(vals[,4]/S.f1, 3)
			norm.f2s <- round(vals[,5]/S.f2, 3)
			norm.f1.gls <- round(vals[,7]/S.f1, 3)
			norm.f2.gls <- round(vals[,8]/S.f2, 3)
			
			if (is.data.frame(return.vals)) {
				new.vals <- data.frame(vals[,1], vals[,2], vals[,3], norm.f1s, norm.f2s, norm.f1.gls, norm.f2.gls)
				names(new.vals) <- c("Speaker", "Vowel", "Context", "F1/S(F1)", "F2/S(F2)", "F1g/S(F1g)", "F2g/S(F2g)")
				return.vals<-rbind(return.vals, new.vals)
			} else {
				return.vals <- data.frame(vals[,1], vals[,2], vals[,3], norm.f1s, norm.f2s, norm.f1.gls, norm.f2.gls)
				names(return.vals) <- c("Speaker", "Vowel", "Context", "F1/S(F1)", "F2/S(F2)", "F1g/S(F1g)", "F2g/S(F2g)")
			}
		}
	}
	attr(return.vals, "norm.method")<-"Watt & Fabricius"
	if (mod.WF) { attr(return.vals, "norm.variant") <- "ModWF" }
	attr(return.vals, "unit.type")<-attr(vowels,"unit.type")
	attr(return.vals, "no.f3s")<-TRUE
	return.vals
}


#
# PLOTTING FUNCTIONS
# 

default.point.colors<- function() {
	c("#ff0000", "#0000ff", "#00ff00", "#ff9900", "#cd00cd", "#00cdcd", "#cdcd00", "#800000", "#008000", "#333399", "#ffff00", "#808080", "#cc99ff")
}

setup.point.color <- function(vowels, color = NA,color.choice=NA){
	if (is.na(color)) {
		pl.c<-"black"	
	} else {
		ind <- 1 # if coloring anything, default to speakers
		if (color == "vowels") ind <- 2
		pl.c<-vector(length=length(vowels[,ind]))
		ind.types <- unique(vowels[,ind])
	
		clrs<-default.point.colors()
		options(warn=-1) # suppressing warnings here 
		if (!is.na(color.choice)) {
			clrs<-color.choice
		}
		options(warn=0)
		clrs<-rep(clrs, len=length(ind.types))
		for (i in 1:length(ind.types)) {
			pl.c[vowels[,ind]==ind.types[i]]<-clrs[i]
		}
	}
	pl.c
}

setup.sizes <- function(vowels, size=NA, a.size=NA, l.size=NA) {
	if (!is.na(size)) {
		p.size<-size
	} else {
		p.size<-0.6+(1/(log(length(vowels[,1]), 10)))
	}
	if (is.na(a.size)) a.size<-1.0
	if (is.na(l.size)) l.size<-a.size*0.75
	c(p.size, a.size, l.size)
}

setup.point.shape <- function(vowels, shape = "speakers", shape.choice = NA) {
	ind <- 1 # default to speakers
	if (shape == "vowels") ind <- 2
	pl.p<-vector(length=length(vowels[,ind]))
	ind.types <- unique(as.character(vowels[,ind]))
	
	#pnts<-c(15, 16, 17, 18, 7, 8, 6, 14)
	pnts<-c(16, 0, 17, 8, 9, 15, 3, 12, 18, 5)
	if (!is.na(shape.choice)) {
		pnts<-shape.choice
	}
	pnts<-rep(pnts, len=length(ind.types))
	for (i in 1:length(ind.types)) {
		pl.p[vowels[,ind]==ind.types[i]]<-pnts[i]
	}
	pl.p
}

setup.axes <- function(vowels) {

  f3.plus <- 0
  if (is.null(attributes(vowels)$no.f3s)) f3.plus <- 1  

  # determine scale for axes
  xstart <- -1000
  xend <- 100000
  ystart <- -1000
  yend <- 100000
  for (i in 1:length(vowels)) {
     if (xstart < max(vowels[,5], vowels[,(7+f3.plus)], na.rm=TRUE)) xstart <- max(vowels[,5], vowels[,(7+f3.plus)], na.rm=TRUE)
     if (xend > min(vowels[,5], vowels[,(7+f3.plus)], na.rm=TRUE)) xend <- min(vowels[,5], vowels[,(7+f3.plus)], na.rm=TRUE)
     if (ystart < max(vowels[,4], vowels[,(6+f3.plus)], na.rm=TRUE)) ystart <- max(vowels[,4], vowels[,(6+f3.plus)], na.rm=TRUE)
     if (yend > min(vowels[,4], vowels[,(6+f3.plus)], na.rm=TRUE)) yend <- min(vowels[,4], vowels[,(6+f3.plus)], na.rm=TRUE)
  }
  xunit <- (xstart-xend)/12
  xstart <- xstart+xunit
  xend <- xend-xunit
  yunit <- (ystart-yend)/12
  ystart <- ystart+yunit
  yend <- yend-yunit
  
  # Reverse axes if the values are in Bark difference
  if (!is.null(attributes(vowels)$norm.method)) {
  	if (attributes(vowels)$norm.method=="Bark Difference") {
  		xend.t<-xend
  		xend<-xstart
  		xstart<-xend.t
  		yend.t<-yend
  		yend<-ystart
  		ystart<-yend.t
  	}
  }

  c(xstart, xend, ystart, yend, xunit, yunit)
}

vowelplot <- function(vowels, speaker = NA, color = NA, color.choice = NA, shape="speakers", shape.choice = NA, size = NA, labels = "none", leg="speakers", a.size = NA, l.size = NA, title = "", subtitle = NA, xlim = NA, ylim = NA) {
  
  ltext<-NA
  nmethod<-"non-"
  if (!is.null(attributes(vowels)$norm.method)) {
	nmethod<-paste(attributes(vowels)$norm.method, " ", sep="")
  }
  vtext<-"Individual"
  if (!is.null(attributes(vowels)$mean.values)) { 
  	vtext<-"Mean"
  } else if (!is.null(attributes(vowels)$median.values)) { 
  	vtext<-"Median"
  } else if (!is.null(attributes(vowels)$standard.devs)) { 
  	vtext<-"Standard deviation of"
  }
  mtext<-paste(vtext, " vowel formant values\n", nmethod, "normalized", sep="")
  
  if (!is.na(speaker)) {
  	vowels<-vowels[vowels[,1]==speaker,]
  	if (labels!="none") ltext<-vowels[,2]  	
  } else if (labels=="vowels") {
  	ltext<-vowels[,2]
  } else if (labels=="speakers") {
  	ltext<-vowels[,1]
  } else if (labels!="none") {
  	if (length(unique(vowels[,1]))==1) {
  		ltext<-vowels[,2]
  	} else {
	  	ltext<-paste(vowels[,1], vowels[,2], sep=":\n")
  	}
  }
  spkrs<-as.character(unique(vowels[,1]))
  
  stext<-""
  if (is.na(subtitle)) {
  	if (!is.null(attributes(vowels)$norm.variant)) {
  		stext<-paste("Variant:", attributes(vowels)$norm.variant)
  	} else if (!is.null(attributes(vowels)$unit.type)) {
  		if (attributes(vowels)$unit.type!="Hertz") {
  			stext<-paste("Units:", attributes(vowels)$unit.type)
  		}
  	}
  } else if (subtitle!="") {
  	stext<-subtitle
  }
  
  axes <- setup.axes(vowels)
  # override computation of axes limits if user specified both
  if (!is.na(xlim[1]) & !is.na(ylim[1])) {
  	axes <- c(xlim[1], xlim[2], ylim[1], ylim[2], (xlim[2]-xlim[1])/12, (ylim[2]-ylim[1])/12)
  }
  pl.c <- setup.point.color(vowels, color, color.choice)
  pl.p <- setup.point.shape(vowels, shape, shape.choice)
  szs <- setup.sizes(vowels, size, a.size, l.size)
  p.s <- szs[1]
  a.s <- szs[2]
  l.s <- szs[3]
  
  if (title != "") mtext <- title
  
  plot(vowels[,5], vowels[,4], xlim=c(axes[1],axes[2]), ylim=c(axes[3],axes[4]), xlab=names(vowels)[5], ylab=names(vowels)[4], pch=pl.p, cex=p.s, cex.main=(a.s + 0.5), cex.axis=(a.s + 0.25), cex.lab=(a.s+0.1), main=mtext, sub=stext, col=pl.c)
  pl.cs<-"black"
  pl.ps<-unique(pl.p)
  if (!(is.null(attributes(vowels)$mean.values) | is.null(attributes(vowels)$median.values))) pl.ps<-unique(pl.p)
  
  if (!is.na(color) & color=="speakers") pl.cs<-unique(pl.c)
  if (!is.na(leg)) {
  	if (leg == "vowels" & (shape == "vowels" | color == "vowels")) {
  	  if (shape=="speakers") pl.ps<-NA
  	  if (!is.na(color) & color=="speakers") {
  	  	pl.lc<-"black"
  	  } else {
  	  	pl.lc<-unique(pl.c)
  	  }
  	  legend("bottomleft", legend=unique(vowels[,2]), col=pl.lc, text.col=pl.lc, pch=pl.ps, inset=.02, cex=l.s, pt.cex=l.s)
  	} else {
  	  if (shape=="vowels") pl.ps<-NA
	  legend("bottomleft", legend=spkrs, col=pl.cs, text.col=pl.cs, pch=pl.ps, inset=.02, cex=l.s)
	}
  }
  f3.plus <- 0
  if (is.null(attributes(vowels)$no.f3s)) {
     f3.plus <- 1
  }
  options(warn=-1) # suppressing warnings here 
  # (arrows create a warning when they're too short to print)
  arrows(vowels[,5], vowels[,4], vowels[,7+f3.plus], vowels[,6+f3.plus], angle=15, length=0.1, lwd=(p.s+0.3), col=pl.c)
  if (!is.na(ltext)) {
	text(vowels[,5], vowels[,4], ltext, adj=c(0,1.5), cex=(2*p.s/3), col=pl.c)
  }
  options(warn=0)
}

add.vowelplot <- function(vowels, speaker=NA, color=NA, color.choice=NA, shape="speakers", shape.choice = NA, size = NA, labels = "none") {
  
  ltext<-NA
  if (!is.na(speaker)) {
  	vowels<-vowels[vowels[,1]==speaker,]
  	if (labels!="none") ltext<-vowels[,2]  	
  } else if (labels=="vowels") {
  	ltext<-vowels[,2]
  } else if (labels=="speakers") {
  	ltext<-vowels[,1]
  } else if (labels!="none") {
  	if (length(unique(vowels[,1]))==1) {
  		ltext<-vowels[,2]
  	} else {
	  	ltext<-paste(vowels[,1], vowels[,2], sep=":\n")
  	}
  }
  spkrs<-as.character(unique(vowels[,1]))
  
  pl.c <- setup.point.color(vowels, color, color.choice)
  pl.p <- setup.point.shape(vowels, shape, shape.choice)
  szs <- setup.sizes(vowels, size)
  p.s <- szs[1]
  a.s <- szs[2]
  l.s <- szs[3]
  
  points(vowels[,5], vowels[,4], pch=pl.p, cex=p.s, cex.lab=(a.s+0.1), main=mtext, col=pl.c)
  pl.cs<-"black"
  pl.ps<-pl.p
  if (!(is.null(attributes(vowels)$mean.values) | is.null(attributes(vowels)$median.values))) pl.ps<-unique(pl.p)
  if (shape=="vowels") pl.ps<-NA
  if (!is.na(color) & color=="speakers") pl.cs<-unique(pl.c)
  f3.plus <- 0
  if (is.null(attributes(vowels)$no.f3s)) {
     f3.plus <- 1
  }
  options(warn=-1) # suppressing warnings here 
  # (arrows create a warning when they're too short to print)
  arrows(vowels[,5], vowels[,4], vowels[,7+f3.plus], vowels[,6+f3.plus], angle=15, length=0.1, lwd=(p.s+0.3), col=pl.c)
  if (!is.na(ltext)) {
	text(vowels[,5], vowels[,4], ltext, adj=c(0,1.5), cex=(2*p.s/3), col=pl.c)
  }
  options(warn=0)
}

add.spread.vowelplot <- function(vowels, mean.points=FALSE, sd.mult=2, ellipsis=FALSE, speaker=NA, color=NA, color.choice=NA, shape="speakers", shape.choice = NA, size = NA, leg=FALSE, labels = "none", separate=TRUE) {
  if (!is.numeric(sd.mult)) { sd.mult<-as.numeric(sd.mult) }

  if (!is.na(speaker)) {
  	vowels<-vowels[vowels[,1]==speaker,]
  }
  spkrs<-as.character(unique(vowels[,1]))
  
  vmns<-compute.means(vowels, separate=separate)
  vsds<-compute.sds(vowels, separate=separate)
  
  ltext<-NA
  if (!is.na(speaker) & (labels != "none")) {
  	ltext<-vmns[,2]  	
  } else if (labels=="vowels") {
  	ltext<-vmns[,2]
  } else if (labels=="speakers") {
  	ltext<-vmns[,1]
  } else if (labels!="none") {
  	if (length(unique(vmns[,1]))==1) {
  		ltext<-vmns[,2]
  	} else {
	  	ltext<-paste(vmns[,1], vmns[,2], sep=":\n")
  	}
  }
  
  pl.c <- setup.point.color(vmns, color, color.choice)
  pl.p <- setup.point.shape(vmns, shape, shape.choice)
  szs <- setup.sizes(vmns, size)
  p.s <- szs[1]+0.5
  a.s <- szs[2]
  l.s <- szs[3]
  
  if (mean.points) points(vmns[,5], vmns[,4], pch=pl.p, cex=p.s, cex.lab=(a.s+0.1), main=mtext, col=pl.c)
  if (ellipsis=="horizvert") {
   	t <- seq(0,7,0.1)
   	# if pl.c is simply "black" need to make it a long enough vec
   	if (length(pl.c) < length(vmns[,5])) {
   		pl.c<-rep(pl.c, length.out=length(vmns[,5]))
   	}
   	for (v in 1:length(vmns[,5])) {
		x <- vmns[v,5] + ((sd.mult*vsds[v,5])*cos(t))
		y <- vmns[v,4] + ((sd.mult*vsds[v,4])*sin(t))
		lines(x, y, lty=2, col=pl.c[v])
	}
  } else if (ellipsis) {
   	t <- seq(0,7,0.1)
  	if (length(pl.c) < length(vmns[,5])) {
   		pl.c<-rep(pl.c, length.out=length(vmns[,5]))
   	}
   	if (separate==TRUE) {
   		for (sp in unique(vmns[,1])) {
   			for (v in as.character(vmns[,2])) {
   				svmns <- vmns[vmns[,1]==sp & vmns[,2]==v,]
   				svsds <- vsds[vsds[,1]==sp & vsds[,2]==v,]
   				svwls <- vowels[vowels[,1]==sp & vowels[,2]==v,]
   				spl.c <- pl.c[vmns[,1]==sp & vmns[,2]==v]
	 			if (nrow(svwls) > 1) {
	   				# Find the principal components and standard devs of the distribution
 		 			prcs <- prcomp(cbind(svwls[,5],svwls[,4]))
	 	 			if (prcs$sdev[1] >= prcs$sdev[2]) {
    					xscale <- prcs$sdev[1] * sd.mult
   					 	yscale <- prcs$sdev[2] * sd.mult
    					column <- 1
  					} else {
   			 			xscale <- prcs$sdev[2] * sd.mult  
    					yscale <- prcs$sdev[1] * sd.mult 
    					column <- 2
  					}
					x <- xscale*cos(t)
					y <- yscale*sin(t)
					rotatedx <- svmns[,5] + x*prcs$rotation[1,column] - y*prcs$rotation[2,column]
					rotatedy <- svmns[,4] + x*prcs$rotation[2,column] + y*prcs$rotation[1,column]
					lines(rotatedx, rotatedy, lty=2, col=spl.c)
				}
			}
		}
	} else {
		for (v in as.character(vmns[,2])) {
			if (!is.na(vsds[vsds[,2]==v,5])) {
	   			# Find the principal components and standard devs of the distribution
 		 		prcs <- prcomp(cbind(vowels[vowels[,2]==v,5],vowels[vowels[,2]==v,4]))
	 	 		if (prcs$sdev[1] >= prcs$sdev[2]) {
    				xscale <- prcs$sdev[1] * sd.mult
   			 		yscale <- prcs$sdev[2] * sd.mult
    				column <- 1
  				} else {
    				xscale <- prcs$sdev[2] * sd.mult  
    				yscale <- prcs$sdev[1] * sd.mult 
    				column <- 2
  				}
				x <- xscale*cos(t)
				y <- yscale*sin(t)
				rotatedx <- vmns[vmns[,2]==v,5] + x*prcs$rotation[1,column] - y*prcs$rotation[2,column]
				rotatedy <- vmns[vmns[,2]==v,4] + x*prcs$rotation[2,column] + y*prcs$rotation[1,column]
				lines(rotatedx, rotatedy, lty=2, col=pl.c[which(vmns[,2]==v)])
			}
		}
	}
  } else {
  	arrows(vmns[,5]-(sd.mult*vsds[,5]), vmns[,4], vmns[,5]+(sd.mult*vsds[,5]), vmns[,4], length=0.1, angle=90, code=3, lty=2, col=pl.c)
 	arrows(vmns[,5], vmns[,4]-(sd.mult*vsds[,4]), vmns[,5], vmns[,4]+(sd.mult*vsds[,4]), length=0.1, angle=90, code=3, lty=2, col=pl.c)
  }
  if (leg) legend("bottomright", legend=paste(sd.mult, " SDs", sep=""), lty=2, cex=0.8, inset=.02)
  options(warn=-1) # suppressing warnings here 
  if (!is.na(ltext)) {
	text(vmns[,5], vmns[,4], ltext, adj=c(0,1.5), cex=(2*p.s/3), col=pl.c)
  }
  options(warn=0)
}