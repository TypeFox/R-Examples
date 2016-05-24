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

#oligoCorr.R

oligoCorr <- function(x, content=c("dna", "rna"), 
case=c("lower", "upper", "as is"), circular=TRUE, disambiguate=TRUE, ...)
{
	if (!is.character(x))
		x <- as.character(x)
#	if (!is.numeric(k) || length(k)!=1 || k<1 || k!=floor(k))#
#		stop("k must be a positive whole number")
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
#Get counts and generate profiles
	profiles <- lapply(2:5, function(k) {
		counts <- switch(as.character(k),
			"2"=pair.counts(x, "as is", circular),
			"3"=triple.counts(x, "as is", circular),
			"4"=quadruple.counts(x, "as is", circular),
			cylinder.counts(x, 1:k, "as is", circular)
		) #switch
		structure(list(content=content, wordLength=k, circular=circular,
		profile=array2vector(counts, sep="", sort=TRUE)), class="OligoProfile")
	}) #lapply
#Generate plots
	par(mfrow=c(2,2))
	for (k in 1:4)
	{
		plot(profiles[[k]], which=3)
	} #for k
	invisible(profiles)
} #function
