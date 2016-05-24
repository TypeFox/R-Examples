#This file is part of the source code for
#SPGS: an R package for identifying statistical patterns in genomic sequences.
#Copyright (C) 2015  Universidad de Chile and INRIA-Chile
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#A copy of Version 2 of the GNU Public License is available in the 
#share/licenses/gpl-2 file in the R installation directory or from 
#http://www.R-project.org/Licenses/GPL-2.

#profile.R

oligoProfile <- function(x, k, content=c("dna", "rna"), 
case=c("lower", "upper", "as is"), circular=TRUE, disambiguate=TRUE, 
plot=TRUE, ...)
{
#Check arguments
	if (!is.character(x))
		x <- as.character(x)
	if (!is.numeric(k) || length(k)!=1 || k<1 || k!=floor(k))
		stop("k must be a positive whole number")
	content <- match.arg(content)
	case <- match.arg(case)
	if (case=="lower") x <- tolower(x)
	if (case=="upper") x <- toupper(x)
	if (!is.logical(circular)) circular <- TRUE
  uniqueSymbols <- as.character(unique(x))
  hasT <- length(intersect(c("t", "T"), uniqueSymbols))>0
  hasU <- length(intersect(c("u", "U"), uniqueSymbols))>0
  if (hasT && hasU)
  {
  	warning(paste(strwrap("the sequence in x contains both thymine \"T\" and uracil \"U\" so correct computation of the complement or reverse complement cannot be assured"), collapse="\n"))
		detectedContent <- content
  }
  else if (hasT && !hasU) detectedContent <- "dna"
  else if (!hasT && hasU) detectedContent <- "rna"
  else detectedContent <- content
  if (content!=detectedContent)
  	warning(paste(strwrap(paste0("since x contains ", detectedContent, " data, the specified content type of \"", content, "\" has been overridden")), collapse="\n"))
	if (disambiguate) x <- disambiguate(x, case="as is")
#Get counts and generate profile
	counts <- switch(as.character(k),
		"2"=pair.counts(x, "as is", circular),
		"3"=triple.counts(x, "as is", circular),
		"4"=quadruple.counts(x, "as is", circular),
		cylinder.counts(x, 1:k, "as is", circular)
	)
	profile <- structure(list(wordLength=k, content=content, case=case, 
	  circular=circular, disambiguate=disambiguate,
	  profile=array2vector(counts, sep="", sort=TRUE)), 
    class="OligoProfile")
	if (!plot) return(profile)
	plot(profile, ...)
	invisible(profile)
} #function

plot.OligoProfile <- function(x, which=1L, units=c("percentage", "count", "proportion"),
main=NULL, xlab=NULL, ylab=NULL, ...)
{
#Check arguments
	if (!is.numeric(which) || which<1 || which>3 || which!=floor(which))
		stop("which must be 1, 2 or 3")
	units <- match.arg(units)
	profile <- switch(units,
		percentage=x$profile/sum(x$profile)*100,
		count=x$profile,
		proportion=x$profile/sum(x$profile)
	)
#Plot profile
	if (which==1L || which==2L)
	{
#Prepare annotation
		if (is.null(main))
		{
			main <- paste0(x$wordLength, "-mer Oligo Profile")
			if (which==1L)
				main=paste(main, "on the Primary Strand")
			else
			{
				main <- paste(main, "on the Complementary Strand")
				profile <- profile[sort(reverseComplement(names(profile), content=x$content))]
			} #if
		} #if
		if (is.null(xlab)) xlab <- paste0(x$wordLength, "-mer")
		if (is.null(ylab))
			ylab <- switch(units, percentage="Percentage", count="Count", proportion="Proportion")
		words <- names(profile)
		if (length(words)>16) pos <- unique(seq.int(1L, length(words), length.out=16L))
		else pos <- seq_along(words)
		words <- words[pos]
		maxWordLen <- max(nchar(words))
		if (maxWordLen>2)
		{
			nDots <- ifelse(maxWordLen<5, maxWordLen-2, 3)
		words <- sapply(words, function(w) {
			paste0(substr(w, 1L, 2L), paste(rep(".", nDots), collapse=""))
		})
		} #if
#Generate plot
		plot(profile, type="l", main=main, xlab=xlab, ylab=ylab, xaxt="n", ...)
		axis(side=1, at=c(0, pos), labels=c("", words), cex=0.6, las=2)
	} #if
	if (which==3L)
	{
#Prepare annotation
		if (is.null(main))
		{
			main <- paste0("Interstrand ", x$wordLength, "-mer Profile Plot")
#			else sub <- ""
		} #if
		if (is.null(xlab))
			xlab <- paste(switch(units, percentage="Percentage", count="Count", proportion="Proportion"), 
				"on Watson Strand")
		if (is.null(ylab)) 
			ylab <- paste(switch(units, percentage="Percentage", count="Count", proportion="Proportion"), 
				"on Crick Strand")
		rcProfile <- profile[reverseComplement(names(profile), content=x$content)]
		lim <- range(profile)
#Generate plot
		plot(profile, rcProfile, type="p", col="red", xlim=lim, ylim=lim, 
#		xaxp=c(lim, 2), yaxp=c(lim, 2), 
			main=main, xlab=xlab, ylab=ylab,
			sub=paste0("Interstrand ", x$wordLength, "-mer correlation:  ", 
			formatC(cor(profile, rcProfile), digits=4L)))
			abline(a=0, b=1, lty=3) #lty=3 is "dotted"
			legend("topleft", legend=c("Corr. 1", "Profiles"), 
			lty=c(3, 0), pch=c(NA, 1), col=c("black", "red"))
	} #if
	invisible()
} #function

print.OligoProfile <- function(x, which=1L, units=c("percentage", "count", "proportion"), 
digits=switch(units, percentage=3L, count=NULL, proportion=3L), ...)
{
#Check arguments
	if (!is.numeric(which) || which<1 || which>2 || which!=floor(which))
		stop("which must be 1 or 2")
	units <- match.arg(units)
	profile <- switch(units,
		percentage=x$profile/sum(x$profile)*100,
		count=x$profile,
		proportion=x$profile/sum(x$profile)
	)
	rcProfile <- profile[reverseComplement(names(profile), content=x$content)]
	ikmercor <- cor(profile, rcProfile)
	main <- paste0(x$wordLength, "-mer Oligo Profile")
	if (which==1L)
		main <- paste(main, "on the Primary Strand")
	else
	{
		main <- paste(main, "on the Complementary Strand")
		profile <- rcProfile[sort(names(rcProfile))]
	} #if
#Print details of NucleotideProfile object
	cat("\n")
	cat(strwrap(main, prefix="\t"), sep="\n")
	cat("\n")
	cat("Content type =", x$content, "\n")
	cat("Circular = ", x$circular, "\n")
	cat("Units =", units, "\n")
	cat("Interstrand", paste0(x$wordLength, "-mer"), "correlation:", ikmercor, "\n")
	print(profile, digits=digits, ...)
	cat("\n")
	invisible(x)
} #function
