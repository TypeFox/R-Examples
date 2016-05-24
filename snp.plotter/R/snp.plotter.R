# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##########################################################################

#' SNP/haplotype association p-value and linkage disequilibrium plotter Creates
#' plots of p-values using single SNP and/or haplotype data. Main features of
#' the package include options to display a linkage disequilibrium (LD) plot and
#' the ability to plot multiple set of results simultaneously. Plots can be
#' created using global and/or individual haplotype p-values along with single
#' SNP p-values. Images are created as either PDF/EPS files.
#' @param config.file Name of a configuration file for snp.plotter parameters in
#'   the form ATTRIBUTE=value. This option can be used in place of specifying 
#'   options at the R command line.
#' @param SNP.FILE Tab-delimited input file containing p-values for single SNPs 
#'   (see note below). The contents of each SNP.FILE includes four necessary
#'   columns ASSOC, SNP.NAME, LOC, and SS.PVAL corresponding to positive or
#'   negative association (indicating susceptibility or protective alleles), a
#'   SNP label, the location, and a p-value for each SNP. SNP labels cannot
#'   start with numbers. Example: "s1.txt,s2.txt,s3.txt" MANDATORY
#' @param HAP.FILE Tab-delimited input file containing p-values for haplotypes 
#'   (see note below). The contents of each HAP.FILE includes three necessary
#'   columns ASSOC, G.PVAL, and I.PVAL followed by a set of columnns of SNPs
#'   with corresponding haplotypes. Haplotypes are presented in a step-wise
#'   fashion with the major allele given as 1 and the minor allele as 2;
#'   haplotype variants for a set of SNPs should be grouped. SNP labels in
#'   HAP.FILE must be the same as in SNP.FILE, and only SNPs with corresponding
#'   haplotypes need to be included. Example: "h1.txt,h2.txt,h3.txt" OPTIONAL
#' @param GENOTYPE.FILE Tab-delimited input file containing genotypes as a PED
#'   file with 6 columns preceding the genotype data: family ID, individual ID,
#'   father ID, mother ID, sex, and affection status; these coloums are not used
#'   in the creation of the LD plot. This file is used for calculating D' or
#'   r-squared values for the LD heatmap plot. Only one LD plot can be shown
#'   (see note below). OPTIONAL
#' @param EVEN.SPACED Logical. Should the p-values be displayed at even spacing
#'   or at genetic map distances?
#' @param USE.GBL.PVAL Logical. Use global haplotype p-values (as opposed to
#'   individual p-values)? Unfilled symbols connected by solid lines are used to
#'   indicate global haplotype p-values, default symbol: circle. Unfilled and
#'   filled symbols are used to indicate alleles 1 and 2, respectively connected
#'   by solid lines and dashed lines for positive and negative association
#'   (indicating susceptibility or protective haplotypes) when using individual
#'   haplotype p-values.
#' @param DISP.HAP Logical. Display haplotype p-values?
#' @param DISP.SNP Logical. Display single SNP p-values?
#' @param DISP.LDMAP Logical. Display the LD heatmap?
#' @param DISP.PHYS.DIST Logical. Display the range of the X-scale?
#' @param DISP.LEGEND Logical. Display a legend with sample labels and
#'   corresponding symbols?
#' @param DISP.COLOR.BAR Logical. Display bar showing colors and corresponding
#'   values of LD plot?
#' @param DISP.TYPE Options: "symbol"
#' @param DISP.MULT.LAB.X Logical. Display evenly spaced X-axis tick-labels; up
#'   to 5 labels are shown.
#' @param DISP.MARKER.LINES Logical. Display lines at p-value thresholds of
#'   0.05, 0.01, 0.001, etc.
#' @param DISP.SNP.NAMES Logical. Display the names of SNPs on a plot.
#' @param DISP.CONNECTING.LINES Logical. Display connecting lines from p-value
#'   plot to LD map.
#' @param USE.COLORS Logical. Restrict LD heatmap colors and default symbol
#'   colors to gray-scale
#' @param COLOR.LIST List of colors (one for each sample) known to GraphApp (see
#'   note below) for displaying p-value symbols. Example:
#'   "red,blue,green,black,orange"
#' @param SYMBOLS Options: circle, square, diamond, triangle; Symbols can either
#'   be filled or not filled by appending "-fill" e.s., square-fill. NA may be
#'   specified. In this case, the SNP.FILE ASSOC column is read and an
#'   up-triangle and down-triangle are used to indicate positive and negative
#'   association (indicating susceptibility or protective alleles),
#'   respectively. Example: "circle,NA,diamond-fill,triangle"
#' @param PALETTE.FILE Colors are hexidecimal HTML color codes; one color per
#'   line. OPTIONAL
#' @param SAMPLE.LABELS Labels for each sample. Example: "d-cc,d2-cc,d1-fam"
#' @param LAB.Y Options: ln (natural log) or log (log10)
#' @param IMAGE.TYPE Options: "pdf" or "eps"
#' @param IMAGE.TITLE Title of the image in quotes. Note: Title text may not
#'   wrap.
#' @param IMAGE.SIZE Options: 3.5 or 7. Sizes are in inches and correspond to 1
#'   or 2 columns per printed page.
#' @param IMAGE.NAME Name of the output file. The correct extension will be
#'   appended depending on the value of IMAGE.TYPE
#' @param PVAL.THRESHOLD The minimum value of the the Y-scale will be set to
#'   this value. Default: 1 (to ignore option).
#' @param LD.TYPE LD metric. Options: "dprime" or "rsquare"
#' @param LD.COLOR.SCHEME LD heatmap color scheme. Options: heat
#'   (red-yellow-white), cm (cyan-magenta), topo (topographical map colors),
#'   gray (gray-scale), or custom; custom requires palette file (PALETTE.FILE)
#'   to be defined
#' @param CONNECTING.LINES.FACTOR Adjusts the length of the connecting lines.
#'   Range: 0-2
#' @param CONNECTING.LINES.ADJ Can be used to adjust the position of connecting
#'   lines in relation to SNP names. Negative values shift the connecting lines
#'   to the left and positive values shift the lines to the right. Range: 0-1
#' @param CONNECTING.LINES.VERT.ADJ Can be used to vertically adjust the
#'   position of connecting lines in relation to SNP names. More negative value
#'   shift the connecting lines down. Range: -0.5-0
#' @param CONNECTING.LINES.FLEX Adjusts the spread of the connecting lines.
#'   Range: 0-2
#' @param SYMBOL.FACTOR Scaling value for symbols. Larger values are generate
#'   larger symbols. Range: 0-1
#' @param FONT.FACTOR Scaling value for SNP names. Larger values are generate
#'   larger SNP names. Range: 0-1
#' @return A list containing two items: config.var and gbl.var, which includes
#'   the values of all significant variables used by snp.plotter
#' @details	snp.plotter produces publishable-quality plots of p-values using
#'   single SNP and/or haplotype data. Main features of the package include
#'   options to display a linkage disequilibrium (LD) plot below the p-value
#'   plot using either the r-squared or D' LD metric with a user-specified LD
#'   heatmap color scheme, setting the X-axis to equal spacing or to use the
#'   physical SNP map, and specification of plot labels, colors and symbols for
#'   desplaying p-values. A major strength of the package is that it can plot
#'   multiple set of results simultaneously. Plots can be created using global
#'   and/or individual haplotype p-values along with single SNP p-values. The
#'   package provides a simple way to convey both association and LD information
#'   in a single appealing graphic and requires virtually no knowledge of the R
#'   programming language. Code to create the LD map was modified from the 
#'   LDHeatmap package by Ji-Hyung Shin, et al. (2006, version 0.2)
#' @note Configuration Files Due to the large number of parameters implemented
#'   for flexibility, it is suggested that snp.plotter be run using the
#'   config.file argument.
#' @note Example Datasets Examples of SNP.FILE, HAP.FILE, GENOTYPE.FILE, and
#'   configuration files are provided at
#'   \url{https://github.com/cannin/snp_plotter} with further
#'   explanation on the file formats.
#' @note Lists Comma delimited lists (SNP.FILE, HAP.FILE, COLOR.LIST, SYMBOLS,
#'   etc) should not have spaces between entries. If using the config.file
#'   argument, these lists should not have quotations in the configuration file.
#'   Example: "red,blue,green,black,orange"
#' @note Colors COLOR.LIST colors are limited to those known to GraphApp. A
#'   short list can be found at
#'   \url{http://en.wikipedia.org/wiki/X11.color.names}; the complete list is
#'   located in the R source code file
#' @note Palettes PALETTE.FILE colors are hexidecimal HTML color codes
#'   \url{http://en.wikipedia.org/wiki/X11.color.names}. The first and last
#'   colors correspond to the lowest and highest value of the chosen LD metric,
#'   respectively. One color per line.
#' @note PDFs The error "unable to start device pdf" may occur when attempting
#'   to overwrite an open PDF document.
#' @note P-values A p-value of 1 or NA can be used in SNP.FILE to prevent
#'   displaying information about a single SNP
#' @note Number of Datasets snp.plotter handles 10 set of results, but provides
#'   default values for only 5 set of results
#' @note File Input All input files should be placed in the same directory
#' @author Augustin Luna \email{augustin.luna at mail.nih.gov}, Kristin K.
#'   Nicodemus \email{kristin.nicodemus at well.ox.ac.uk}. Website:
#'   \url{https://github.com/cannin/snp_plotter}
#' @keywords aplot, hplot
#' @examples
#' \dontrun{
#' snp.plotter(config.file="config.txt")
#' }
#' @import genetics
#' @export
snp.plotter <- function(EVEN.SPACED = FALSE, 
		PVAL.THRESHOLD = 1, 
		USE.GBL.PVAL = TRUE, 
		SYMBOLS = NA, 
		SAMPLE.LABELS = NULL, 
		LAB.Y = "log", 
		DISP.HAP = FALSE, 
		DISP.SNP = TRUE, 
		DISP.COLOR.BAR = TRUE, 
		DISP.PHYS.DIST = TRUE, 
		DISP.LEGEND = TRUE, 
		DISP.MARKER.LINES = TRUE, 
		DISP.LDMAP = FALSE, 
		DISP.TYPE = "symbol", 
		DISP.MULT.LAB.X = FALSE, 
		DISP.SNP.NAMES = TRUE,
		DISP.CONNECTING.LINES = TRUE,
		LD.TYPE = "dprime", 
		LD.COLOR.SCHEME = "heat", 
		USE.COLORS = TRUE, 
		COLOR.LIST = NULL, 
		PALETTE.FILE = NULL, 
		IMAGE.TITLE = NULL, 
		IMAGE.NAME = "snp.plotter", 
		IMAGE.TYPE = "pdf", 
		IMAGE.SIZE = 3.5,  
		CONNECTING.LINES.FACTOR = 1, 
		CONNECTING.LINES.ADJ = 0, 
		CONNECTING.LINES.VERT.ADJ = -1, 
		CONNECTING.LINES.FLEX = 0, 
		SNP.FILE = NULL, 
		HAP.FILE = NULL, 
		GENOTYPE.FILE = NULL,
		FONT.FACTOR = NULL,
		SYMBOL.FACTOR = NULL,
		config.file = NULL) {
	
	#Read in configuration file
	read.config <- function(config.file, config.var) {
	
		#DEBUG STATEMENT
		cat("START READ.CONFIG\n")
	
		raw.config <- read.table(config.file, header=FALSE, as.is=TRUE, fill=TRUE, sep="\t", blank.lines.skip = TRUE)
	
		for(i in 1:length(raw.config$V1)) {
			tmp <- strsplit(raw.config$V1[i], "=")
	
			if(identical("SNP.FILE", tmp[[1]][1])) {
				config.var$SNP.FILE <- tmp[[1]][2]
			} 
		
			if(identical("HAP.FILE", tmp[[1]][1])) {
				config.var$HAP.FILE <- tmp[[1]][2]
			} 
	
			if(identical("PALETTE.FILE", tmp[[1]][1])) {
				config.var$PALETTE.FILE <- tmp[[1]][2]
			}
	
			if(identical("LAB.Y", tmp[[1]][1])) {
				config.var$LAB.Y <- tmp[[1]][2]
			} 
	
			if(identical("SYMBOLS", tmp[[1]][1])) {
				config.var$SYMBOLS <- tmp[[1]][2]
			} 
	
			if(identical("EVEN.SPACED", tmp[[1]][1])) {
				config.var$EVEN.SPACED <- as.logical(tmp[[1]][2])
			} 
	
			if(identical("PVAL.THRESHOLD", tmp[[1]][1])) {
				config.var$PVAL.THRESHOLD <- as.numeric(tmp[[1]][2])
			} 
	
			if(identical("USE.GBL.PVAL", tmp[[1]][1])) {
				config.var$USE.GBL.PVAL <- as.logical(tmp[[1]][2])
			} 
	
			if(identical("DISP.TYPE", tmp[[1]][1])) {
				config.var$DISP.TYPE <- tmp[[1]][2]
			} 
		
			if(identical("DISP.LDMAP", tmp[[1]][1])) {
				config.var$DISP.LDMAP <- as.logical(tmp[[1]][2])
			} 
	
			if(identical("DISP.HAP", tmp[[1]][1])) {
				config.var$DISP.HAP <- as.logical(tmp[[1]][2])
			}
			
			if(identical("DISP.SNP", tmp[[1]][1])) {
				config.var$DISP.SNP <- as.logical(tmp[[1]][2])
			}
	
			if(identical("LD.COLOR.SCHEME", tmp[[1]][1])) {
				config.var$LD.COLOR.SCHEME <- tmp[[1]][2]
			} 
	
	        if(identical("USE.COLORS", tmp[[1]][1])) {
				config.var$USE.COLORS <- as.logical(tmp[[1]][2])
			}
	
			if(identical("COLOR.LIST", tmp[[1]][1])) {
				config.var$COLOR.LIST <- tmp[[1]][2]
			} 
	
			if(identical("LD.TYPE", tmp[[1]][1])) {
				config.var$LD.TYPE <- tmp[[1]][2]
			} 
	
			if(identical("DISP.COLOR.BAR", tmp[[1]][1])) {
				config.var$DISP.COLOR.BAR <- as.logical(tmp[[1]][2])
			} 

			if(identical("DISP.PHYS.DIST", tmp[[1]][1])) {
				config.var$DISP.PHYS.DIST <- as.logical(tmp[[1]][2])
			} 

			if(identical("DISP.SNP.NAMES", tmp[[1]][1])) {
				config.var$DISP.SNP.NAMES <- as.logical(tmp[[1]][2])
			} 
	
			if(identical("DISP.CONNECTING.LINES", tmp[[1]][1])) {
				config.var$DISP.CONNECTING.LINES <- as.logical(tmp[[1]][2])
			} 	
	
			if(identical("GENOTYPE.FILE", tmp[[1]][1])) {
				config.var$GENOTYPE.FILE <- tmp[[1]][2]
			} 
	
			if(identical("IMAGE.TITLE", tmp[[1]][1])) {
				config.var$IMAGE.TITLE <- tmp[[1]][2]
			} 
	
			if(identical("DISP.LEGEND", tmp[[1]][1])) {
				config.var$DISP.LEGEND <- as.logical(tmp[[1]][2])
			} 
			
			if(identical("SAMPLE.LABELS", tmp[[1]][1])) {
				config.var$SAMPLE.LABELS <- tmp[[1]][2]
			} 
	
			if(identical("IMAGE.TYPE", tmp[[1]][1])) {
				config.var$IMAGE.TYPE <- tmp[[1]][2]
			}

			if(identical("IMAGE.SIZE", tmp[[1]][1])) {
				config.var$IMAGE.SIZE <- tmp[[1]][2]
			}
	    
			if(identical("DISP.MULT.LAB.X", tmp[[1]][1])) {
				config.var$DISP.MULT.LAB.X <- as.logical(tmp[[1]][2])
			}

			if(identical("IMAGE.NAME", tmp[[1]][1])) {
				config.var$IMAGE.NAME <- tmp[[1]][2]
			}
	
			if(identical("CONNECTING.LINES.FACTOR", tmp[[1]][1])) {
				config.var$CONNECTING.LINES.FACTOR <- as.numeric(tmp[[1]][2])
			}
	
			if(identical("CONNECTING.LINES.ADJ", tmp[[1]][1])) {
				config.var$CONNECTING.LINES.ADJ <- as.numeric(tmp[[1]][2])
			}
			
			if(identical("CONNECTING.LINES.VERT.ADJ", tmp[[1]][1])) {
				config.var$CONNECTING.LINES.VERT.ADJ <- as.numeric(tmp[[1]][2])
			}
	
			if(identical("CONNECTING.LINES.FLEX", tmp[[1]][1])) {
				config.var$CONNECTING.LINES.FLEX <- as.numeric(tmp[[1]][2])
			}
			
			if(identical("FONT.FACTOR", tmp[[1]][1])) {
				config.var$FONT.FACTOR <- as.numeric(tmp[[1]][2])
			}
			
			if(identical("SYMBOL.FACTOR", tmp[[1]][1])) {
				config.var$SYMBOL.FACTOR <- as.numeric(tmp[[1]][2])
			}
		}
	
		if(is.null(config.var$SNP.FILE)) {
			stop("Invalid SNP data file: ", config.var$SNP.FILE, "\n")
		}
	
		if(is.null(config.var$HAP.FILE) & config.var$DISP.HAP) {
			stop("Invalid haplotype data file: ", config.var$HAP.FILE, "\n")
		}
	
		if(is.null(config.var$GENOTYPE.FILE) & config.var$DISP.LDMAP) {
			stop("Invalid genotype data file: ", config.var$GENOTYPE.FILE, "\n")
		}
	
		#DEBUG STATEMENT
		cat("FINISH READ.CONFIG\n")
	
		return(config.var)
	}
	
	retrieve.data <- function(split.snp.file, split.hap.file, config.var, gbl.var) {
	
		#DEBUG STATEMENT
		cat("START RETRIEVE.DATA\n")

		#initialize distance variables
		min.dist <- NULL
		max.dist <- NULL
		
		for(i in 1:length(split.snp.file[[1]])) {

			gbl.var$snp.data <- read.delim(split.snp.file[[1]][i], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)

			gbl.var$snp.data$SS.PVAL[which(is.na(gbl.var$snp.data$SS.PVAL))] <- 1

			#DEBUG STATEMENT
			#cat("gbl.var$snp.data ", length(gbl.var$snp.data), "\n")
	
			is.complete.snp.file <- TRUE
	
			#check for data columns and compatibility with display type option
			if(is.null(gbl.var$snp.data$MAJOR.ALLELE) & ((config.var$DISP.TYPE == "allele") | (config.var$DISP.TYPE == "base"))) {
				stop("Missing SNP data file column MAJOR.ALLELE. Neither DISP.TYPE \"allele\" nor \"base\" can be chosen.\n")
			} else if (is.null(gbl.var$snp.data$BASE) & ((config.var$DISP.TYPE == "allele") | (config.var$DISP.TYPE == "base"))) {
				stop("Missing SNP data file column BASE. Neither DISP.TYPE \"allele\" nor \"base\" can be chosen.\n")
			} else if (is.null(gbl.var$snp.data$ASSOC) & ((config.var$DISP.TYPE == "allele") | (config.var$DISP.TYPE == "base"))) {
				stop("Missing SNP data file column ASSOC. Neither DISP.TYPE \"allele\" nor \"base\" can be chosen.\n")
			} else if (is.null(gbl.var$snp.data$ASSOC) & ((config.var$DISP.TYPE == "symbol") & (any(is.na(gbl.var$symbol.list))))) {
				stop("Missing SNP data file column ASSOC. This column is necessary when symbol type is specified as NA\n")
			} else if (is.null(gbl.var$snp.data$SNP.NAME)) {
				stop("Missing SNP data file column SNP.NAME\n")
			} else if (is.null(gbl.var$snp.data$LOC)) {
				stop("Missing SNP data file column LOC\n")
			} else if (is.null(gbl.var$snp.data$SS.PVAL)) {
				stop("Missing SNP data file column SS.PVAL\n")
			} else if (length(grep("^[0-9]", gbl.var$snp.data$SNP.NAME)) > 0) {
				stop("SNP names cannot start with numbers.\n")
			} else if (is.null(gbl.var$snp.data$ASSOC) & length(grep("NA", config.var$SYMBOLS)) > 0) {
				stop("Missing SNP data file column ASSOC. The column is required when using symbol type 'NA'\n")
			}
	
			#determine is a simplified sample set has been inputted
			if(is.null(gbl.var$snp.data$MAJOR.ALLELE) | is.null(gbl.var$snp.data$BASE) | is.null(gbl.var$snp.data$ASSOC)) {
	            		is.complete.snp.file <- FALSE
			}
			
			if(config.var$DISP.HAP) {
				gbl.var$hap.data <- read.delim(split.hap.file[[1]][i], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE)
			}
	
			cur.snp.num <- nrow(gbl.var$snp.data)
	
			#create 2 hash tables containing SNP names and positions for current and previous sample sets
			#one hash has the name as the key and the other the location because R has no way of looking
			#up hash keys based on values
			for(p in 1:cur.snp.num) {
	
				if(!is.na(gbl.var$snp.data$SNP.NAME[p]) | !is.na(gbl.var$snp.data$LOC[p])) {
					assign(gbl.var$snp.data$SNP.NAME[p], list(position=gbl.var$snp.data$LOC[p]), envir=gbl.var$snp.hash.names.pos)
	    	        assign(as.character(gbl.var$snp.data$LOC[p]), list(snp.name=gbl.var$snp.data$SNP.NAME[p]), envir=gbl.var$snp.hash.pos.names)
				}
			}
	
			#create a variable with only the SNP names
	        snp.hash.names <- ls(envir=gbl.var$snp.hash.names.pos)
			#create a variable with only the SNP positions
			snp.hash.pos <- ls(envir=gbl.var$snp.hash.pos.names)
	
			#DEBUG STATEMENT 
			#cat("length(snp.hash.names) ", length(snp.hash.names), "\n")
			#cat("length(snp.hash.pos) ", length(snp.hash.pos), "\n")

	        #create a variable with only the sorted SNP positions
			gbl.var$sorted.snp.pos <- sort(as.numeric(snp.hash.pos))
	
			#sort the SNP names
	        for(n in 1:length(snp.hash.pos)) {
				#cat("POS: ", snp.hash.pos[i], "NAME: ", get(snp.hash.pos[i], env=snp.hash.pos.name)$name, "\n")
				gbl.var$sorted.snp.names[n] <- get(snp.hash.pos[n], envir=gbl.var$snp.hash.pos.names)$snp.name
			}
	
			#get the number of unique SNP in all current sample sets 
			gbl.var$snp.num <- length(gbl.var$sorted.snp.names)

			#DEBUG STATEMENT
			#cat("gbl.var$snp.num ", gbl.var$snp.num, "\n")

			#Set the gbl.var$cex.factor, which determines font size
			formatting.var.list <- set.image.parameters(config.var, gbl.var)

			gbl.var$font.size <- formatting.var.list$font.size
			gbl.var$line.width <- formatting.var.list$line.width
			gbl.var$cex.factor <- formatting.var.list$cex.factor
			gbl.var$cex.factor.symbol <- formatting.var.list$cex.factor.symbol

			#DEBUG STATEMENT
			#cat("gbl.var$cex.factor ", gbl.var$cex.factor, "\n")
			#cat("gbl.var$snp.data$BASE ", nrow(gbl.var$snp.data$BASE), "\n")
		
			if(is.complete.snp.file) {
				gbl.var$snp.bases <- strsplit(gbl.var$snp.data$BASE, "/")
			}
	
			if(config.var$IMAGE.SIZE == 3.5) {
				gbl.var$snp.names <- substr(gbl.var$snp.data$SNP.NAME, 1, 9)
			} else if(config.var$IMAGE.SIZE == 7) {
				gbl.var$snp.names <- substr(gbl.var$snp.data$SNP.NAME, 1, 9)
			} else {
				stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
			}
	
			if(config.var$DISP.LDMAP & (i == 1)) {
	
				if(!is.null(config.var$GENOTYPE.FILE)) {
					geno.data.raw <- read.table(config.var$GENOTYPE.FILE, header=FALSE, blank.lines.skip = TRUE)
				} else {
					stop("Missing GENOTYPE.FILE from configuration file\n")
				}	

				#DEBUG STATEMENT
				#cat("geno.data.raw ", ncol(geno.data.raw), "\n")

				snp.headers <- c("fam", "ind", "fid", "mid", "sex", "aff")
	
				for(i in 1:nrow(gbl.var$snp.data)) {
					snp.headers <- c(snp.headers, paste(gbl.var$snp.data$SNP.NAME[i], ".1", sep=""), paste(gbl.var$snp.data$SNP.NAME[i], ".2", sep=""))
				}
	
				for(k in seq(7, length(snp.headers), by=2)) {
					gbl.var$geno.data <- cbind(gbl.var$geno.data, cbind(paste(geno.data.raw[,k], geno.data.raw[,k+1], sep="/")))				

					#DEBUG STATEMENT
					#cat("geno.data", ncol(gbl.var$geno.data), "\n")
				}
	
				for(m in 1:gbl.var$snp.num) {
					gbl.var$geno.data[which(gbl.var$geno.data[,m] == "0/0"), m] <- NA
	
					#DEBUG STATEMENT
					#cat("geno.data", ncol(gbl.var$geno.data), "\n")				
				}
			}
				
			min.dist <- min(c(min.dist, gbl.var$snp.data$LOC), na.rm = TRUE)
			max.dist <- max(c(max.dist, gbl.var$snp.data$LOC), na.rm = TRUE)
			gbl.var$total.dist <- max.dist - min.dist 
	
			fix.var <- fix.values(config.var, gbl.var)
			
			config.var <- fix.var$config.var
			gbl.var <- fix.var$gbl.var
	
			if(config.var$EVEN.SPACED) {
				gbl.var$min.x <- 1
				gbl.var$max.x <- length(gbl.var$sorted.snp.names)
			}
			else {
				gbl.var$min.x <- min(c(gbl.var$min.x, gbl.var$snp.data$LOC), na.rm = TRUE)
				#cat("gbl.var$min.x", gbl.var$min.x, "\n")
				gbl.var$max.x <- max(c(gbl.var$max.x, gbl.var$snp.data$LOC), na.rm = TRUE)
				#cat("gbl.var$max.x", gbl.var$max.x, "\n")
			}
	
			if(config.var$DISP.HAP) {
				if(config.var$USE.GBL.PVAL) {
					gbl.var$min.y <- min(c(gbl.var$min.y, min(gbl.var$hap.data$G.PVAL, na.rm = TRUE), min(gbl.var$snp.data$SS.PVAL, na.rm = TRUE)))
					# Ensure that the y scale extends beyond the data region, round to a whole number
					gbl.var$max.y <- ceiling(max(c(gbl.var$max.y, max(gbl.var$hap.data$G.PVAL, na.rm = TRUE), max(gbl.var$snp.data$SS.PVAL, na.rm = TRUE))))
				}
	
				if(!config.var$USE.GBL.PVAL) {
					gbl.var$min.y <- min(c(gbl.var$min.y, min(gbl.var$hap.data$I.PVAL, na.rm = TRUE), min(gbl.var$snp.data$SS.PVAL, na.rm = TRUE)))
					# Ensure that the y scale extends beyond the data region, round to a whole number
					gbl.var$max.y <- ceiling(max(c(gbl.var$max.y, max(gbl.var$hap.data$I.PVAL, na.rm = TRUE), max(gbl.var$snp.data$SS.PVAL, na.rm = TRUE))))

					#DEBUG STATEMENT
					#cat("gbl.var$hap.data$I.PVAL ", (gbl.var$hap.data$I.PVAL), "\n")
				}
			}
			else {
				gbl.var$min.y <- min(c(gbl.var$min.y, min(gbl.var$snp.data$SS.PVAL, na.rm = TRUE)))
				# Ensure that the y scale extends beyond the data region, round to a whole number
				gbl.var$max.y <- ceiling(max(c(gbl.var$max.y, max(gbl.var$snp.data$SS.PVAL, na.rm = TRUE))))
			}

			#DEBUG STATEMENT
			#cat("gbl.var$max.y ", gbl.var$max.y, "\n")
	
			gbl.var$cur.exp <- 0
	
			max.min.diff <- gbl.var$max.x - gbl.var$min.x
			exp.test <- TRUE
	
			while(exp.test) {
				if(max.min.diff > 10) {
					gbl.var$cur.exp <- gbl.var$cur.exp + 1
					max.min.diff <- max.min.diff / 10
				}
				else {
					exp.test <- FALSE
				}
			}
	
			if(config.var$EVEN.SPACED) {
				gbl.var$r.x <- c((gbl.var$min.x - 0.5), gbl.var$max.x) # 0 - max(LOC)
				#cat(gbl.var$r.x)
			}
			else {
				gbl.var$r.x <- c(gbl.var$min.x - 100, gbl.var$max.x) # 0 - max(LOC)
				#cat(gbl.var$r.x, "\n")
			}
	
			if(config.var$EVEN.SPACED) {
				gbl.var$axis.x <- seq(gbl.var$min.x, gbl.var$max.x)
			}
			else {
				gbl.var$axis.x <- seq(gbl.var$min.x, gbl.var$max.x, 10^(gbl.var$cur.exp - 1)) # BY Variable
				gbl.var$axis.x[length(gbl.var$axis.x)] <- gbl.var$max.x
			}
	
			#DEBUG STATEMENT
			#cat("config.var$DISP.LEGEND ", config.var$DISP.LEGEND, "\n")
			#cat("config.var$DISP.LDMAP ", config.var$DISP.LDMAP, "\n")
			#cat("gbl.var$max.y ", gbl.var$max.y, "\n")

			if(config.var$DISP.LEGEND & !config.var$DISP.LDMAP) {

				gbl.var$axis.y <- seq(0, gbl.var$max.y, 1) # BY Variable
				gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y + 0.25) # 0 - max(PVAL) for either SS or HAP
			}
			else {
				gbl.var$axis.y <- seq(0, gbl.var$max.y, 1) # BY Variable
				gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y) # 0 - max(PVAL) for either SS or HAP
			}
		}

		gbl.var <- draw.plot.grid.setup(config.var, gbl.var)	
		
		for(i in 1:gbl.var$snp.samples) {
	
			gbl.var$cur.sample <- i 

			gbl.var$snp.data <- read.delim(split.snp.file[[1]][i], header=TRUE, sep="\t", as.is=TRUE)
		
			if(config.var$DISP.HAP) {
				gbl.var$hap.data <- read.delim(split.hap.file[[1]][i], header=TRUE, sep="\t", as.is=TRUE)
			}
	
			if(config.var$EVEN.SPACED) {
				for(j in 1:length(gbl.var$sorted.snp.names)) {	
					for(k in 1:length(gbl.var$snp.data$LOC)) {
						if(gbl.var$snp.data$LOC[k] == gbl.var$sorted.snp.pos[j]) {
							gbl.var$snp.data$LOC[k] <- j
						}
					}
				}
			}
			
			fix.var <- fix.values(config.var, gbl.var)
			
			config.var <- fix.var$config.var
			gbl.var <- fix.var$gbl.var
	
			if(gbl.var$cur.sample == 1) {
	
				#DEBUG STATEMENT
				#cat("gbl.var$cur.sample", gbl.var$cur.sample, "\n")
				#cat("length(gbl.var$sorted.snp.names) ", length(gbl.var$sorted.snp.names), "\n")
				#cat("length(gbl.var$sorted.snp.pos) ", length(gbl.var$sorted.snp.pos), "\n")
				#cat("DISP.LDMAP ", config.var$DISP.LDMAP, "\n")
		
				draw.plot.grid.snp.names(config.var, gbl.var)
			}

			if(config.var$DISP.SNP) {
				draw.plot.grid.ss(config.var, gbl.var)
			}
			
			#DEBUG STATEMENT
			#cat("gbl.var$hap.samples ", gbl.var$hap.samples, "\n")
			#cat("config.var$DISP.HAP ", config.var$DISP.SNP, "\n")
	
			if(config.var$DISP.HAP & (gbl.var$hap.samples >= i)) {
				draw.plot.grid.hap(config.var, gbl.var)
			}
		}
	
		#DEBUG STATEMENT
		cat("FINISH RETRIEVE.DATA\n")

		return(list(config.var = config.var, gbl.var = gbl.var))
	}
	
	fix.values <- function(config.var, gbl.var) {
	
		#DEBUG STATEMENT
		cat("START FIX.VALUES\n")
	
		if(config.var$LAB.Y == "ln") {
	
			gbl.var$snp.data$SS.PVAL <- -log(gbl.var$snp.data$SS.PVAL) # LOG Variable
	
			if(config.var$DISP.HAP) {
				gbl.var$hap.data$G.PVAL <- -log(gbl.var$hap.data$G.PVAL) # LOG Variable
				gbl.var$hap.data$I.PVAL <- -log(gbl.var$hap.data$I.PVAL) # LOG Variable
			}
	
	
			#DEBUG STATEMENT 
			#cat("PVAL.THRESHOLD ", config.var$PVAL.THRESHOLD, "\n")
		
			if(gbl.var$pval.flag) {
				if(config.var$PVAL.THRESHOLD != 1) {
					config.var$PVAL.THRESHOLD <- -log(config.var$PVAL.THRESHOLD)
					gbl.var$pval.flag <- FALSE
				} else {
					config.var$PVAL.THRESHOLD <- -1
					gbl.var$pval.flag <- FALSE
				}
			}
		} else if (config.var$LAB.Y == "log") {
			
			gbl.var$snp.data$SS.PVAL <- -log10(gbl.var$snp.data$SS.PVAL) # LOG Variable
	
			if(config.var$DISP.HAP) {
				gbl.var$hap.data$G.PVAL <- -log10(gbl.var$hap.data$G.PVAL) # LOG Variable
				gbl.var$hap.data$I.PVAL <- -log10(gbl.var$hap.data$I.PVAL) # LOG Variable
			}
	
	
			#DEBUG STATEMENT 
			#cat("PVAL.THRESHOLD ", config.var$PVAL.THRESHOLD, "\n")
		
			if(gbl.var$pval.flag) {
				if(config.var$PVAL.THRESHOLD != 1) {
					config.var$PVAL.THRESHOLD <- -log10(config.var$PVAL.THRESHOLD)
					gbl.var$pval.flag <- FALSE
				} else {
					config.var$PVAL.THRESHOLD <- -1
					gbl.var$pval.flag <- FALSE
				}
			}
		} else {
			stop("Invalid LAB.Y: ", config.var$LAB.Y, ". Specify either 'log' or 'ln'.\n")
		}
	
	
		if(config.var$DISP.LDMAP) {
	
			gbl.var$geno.data <- data.frame(gbl.var$geno.data)
	
			for(i in 1:ncol(gbl.var$geno.data)) {
			
				gbl.var$geno.data[,i] <- genotype(gbl.var$geno.data[,i])
			}
		}
	
		#DEBUG STATEMENT
		cat("FINISH FIX.VALUES\n")
	
		return(list(config.var = config.var, gbl.var = gbl.var))
	}
	
	draw.plot.grid.setup <- function(config.var, gbl.var) { 
	
		#DEBUG STATEMENT
		cat("START DRAW.PLOT.GRID.SETUP\n")
	
		#set y-axis initial value	
		if(config.var$PVAL.THRESHOLD != -1) {
			gbl.var$r.y[1] <- config.var$PVAL.THRESHOLD
			gbl.var$axis.y <- seq(floor(gbl.var$r.y[1]), max(gbl.var$axis.y))
		} else {
			gbl.var$r.y[1] <- 0
		}

		#-------------------VIEWPORT BEGINS------------------
		#-------------------TOP VIEWPORT---------------------	
		if(config.var$EVEN.SPACED) {
			top.vp.ldmap <- viewport(layout = grid.layout(5, 3, widths = c(0.17, 0.6875, 0.1425), 
									 heights = c(.05, .4, .1, .4, .05)), 
									 name = "top.vp.ldmap")
		} else {
			top.vp.ldmap <- viewport(layout = grid.layout(5, 3, widths = c(0.18, 0.6875, 0.1325), 
									 heights = c(.05, .4, .1, .4, .05)), 
									 name = "top.vp.ldmap")
		}

		top.vp.noldmap <- viewport(width = 0.675, layout = grid.layout(3, 1, widths = 1, 
	                               heights = c(.1, .575, .225)), 
								   name = "top.vp.noldmap")

		#-------------------DATA VIEWPORT---------------------
		data.vp.ldmap <- viewport(xscale = c(gbl.var$r.x[1], gbl.var$axis.x[length(gbl.var$axis.x)]), 
								  yscale = c(gbl.var$axis.y[1], gbl.var$r.y[2]), 
								  height = 0.8, 
								  width = 0.8, 
								  layout.pos.row = 2, 
								  layout.pos.col = 2, 
								  name = "data.vp.ldmap")
		data.vp.noldmap <- viewport(xscale = c(gbl.var$r.x[1], gbl.var$axis.x[length(gbl.var$axis.x)]), 
									yscale = c(gbl.var$axis.y[1], gbl.var$r.y[2]), 
									height = .8, 
									width = .8, 
									layout.pos.row = 2, 
									layout.pos.col = 1, 
									name = "data.vp.noldmap")
	
		#-------------------TITLE VIEWPORT---------------------
		title.vp.ldmap <- viewport(height = .8, 
								   width = .8, 
								   layout.pos.row = 1, 
								   layout.pos.col = 2, 
								   name = "title.vp.ldmap")
		title.vp.noldmap <- viewport(height = .8, 
									 width = .8, 
									 layout.pos.row = 1, 
									 layout.pos.col = 1, 
								     name = "title.vp.noldmap")
	
		if(config.var$DISP.LDMAP) {
			pushViewport(vpTree(top.vp.ldmap, vpList(data.vp.ldmap)))
		}
		else {
			pushViewport(vpTree(top.vp.noldmap, vpList(data.vp.noldmap)))
		}
	
		#-------------------AXIS BEGINS------------------
	
		image.title.text <- textGrob(config.var$IMAGE.TITLE, x = 0.5, y = 1.2, just = c("center"), 
									 default.units = "native", gp = gpar(fontsize = gbl.var$font.size + 4, fontface = "bold"))
		image.title <- gTree(children=gList(image.title.text), vp=title.vp.ldmap, name="image.title") 
	
		grid.draw(image.title)

		if(config.var$LAB.Y == "ln") {
			if(config.var$DISP.MARKER.LINES) {
	
				if(-log(0.05) > gbl.var$axis.y[1]) {
					grid.lines(x = gbl.var$axis.x, y = -log(0.05), default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width))
				}
				
				initial.abline <- 0.1 

				while(-log(initial.abline) < gbl.var$axis.y[1]) {
					initial.abline <- initial.abline / 10 
					
					#DEBUG STATEMENT
					#cat("initial.abline 1", initial.abline, "\n")
				}
				
				while(-log(initial.abline) < max(gbl.var$axis.y)) {

					#DEBUG STATEMENT
					#cat("initial.abline 2", initial.abline, "\n")

					grid.lines(x = gbl.var$axis.x, y = -log(initial.abline), default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width))
					initial.abline <- initial.abline / 10
				}
			}
		} else {
			if(config.var$DISP.MARKER.LINES) {
				
				if(-log10(0.05) > gbl.var$axis.y[1]) {
					grid.lines(x = gbl.var$axis.x, y = -log10(0.05), default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width))
				}
				
				initial.abline <- 0.1

				while(-log10(initial.abline) < gbl.var$axis.y[1]) {
					initial.abline <- initial.abline / 10 
				}
				
				while(-log10(initial.abline) < max(gbl.var$axis.y)) {
					grid.lines(x = gbl.var$axis.x, y = -log10(initial.abline), default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width))
					initial.abline <- initial.abline / 10
				}
			}
		}
	
		if(config.var$DISP.MULT.LAB.X) {
			#Truncate the number of X-axis labels
			if(length(gbl.var$axis.x) > 5) {
				increment <- ceiling(length(gbl.var$axis.x) / 5)
	  
				axis.x.labels <- gbl.var$axis.x[seq(1, length(gbl.var$axis.x), by= increment)]
	  
				if(length(gbl.var$axis.x) %% 5 != 0) {
					axis.x.labels <- c(axis.x.labels, gbl.var$axis.x[length(gbl.var$axis.x)])
				}
			}
		}
		else {
			#Truncate the number of X-axis labels
			if(length(gbl.var$axis.x) >= 2) {
	  	
				axis.x.labels <- c(head(gbl.var$axis.x, 1), tail(gbl.var$axis.x, 1))
			}
		}

		if(config.var$DISP.LDMAP) {
			if(config.var$EVEN.SPACED) {
				grid.xaxis(at = gbl.var$axis.x, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
			}
			else {
				grid.xaxis(at = axis.x.labels, label = FALSE, gp = gpar(cex = gbl.var$cex.factor), name = "axis.x.labels")
				grid.xaxis(at = gbl.var$sorted.snp.pos, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
				grid.text(axis.x.labels[1], x = -0.025, y = -0.1, just = "right", gp = gpar(fontsize = gbl.var$font.size))
				grid.text(axis.x.labels[2], x = 1.015, y = -0.1, just = "left", gp = gpar(fontsize = gbl.var$font.size))
	
				#Remove trailing tick mark
				grid.edit("axis.x.labels", gp = gpar(col = "white"))
			}
		}
		else {
			if(config.var$EVEN.SPACED) {
				grid.xaxis(at = gbl.var$axis.x, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
			} else {
				grid.xaxis(at = axis.x.labels, label = FALSE, gp = gpar(cex = gbl.var$cex.factor), name = "axis.x.labels")
				grid.xaxis(at = gbl.var$sorted.snp.pos, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
				grid.text(axis.x.labels[1], x = -0.025, y = -0.09, just = "right", gp = gpar(fontsize = gbl.var$font.size))
				grid.text(axis.x.labels[2], x = 1.015, y = -0.09, just = "left", gp = gpar(fontsize = gbl.var$font.size))
	
				#Remove trailing tick mark
				grid.edit("axis.x.labels", gp = gpar(col = "white"))
			}
		}	
	
		grid.yaxis(at = gbl.var$axis.y, label = gbl.var$axis.y, gp = gpar(fontsize = gbl.var$font.size, lwd = gbl.var$line.width))

		#display the physical distance on plots without LD maps 
		if(config.var$DISP.PHYS.DIST & !config.var$DISP.LDMAP) {
		
			if(gbl.var$total.dist > 1000) {
				grid.text(paste("Physical Distance: ", round(gbl.var$total.dist/1000, 1), " kb", sep=""),
						  x = 0.5, 
                          y = -0.2625, 
                          just = "center", 
                          gp = gpar(fontsize = gbl.var$font.size))
			} else {
				grid.text(paste("Physical Distance: ", 
                          gbl.var$total.dist, " bases", sep=""),
						  x = 0.5, y = -0.2625, 
                          just = "center", 
                          gp = gpar(fontsize = gbl.var$font.size))
			}
		}
	
		if(config.var$LAB.Y == "ln") {
			grid.text("-ln(p-value)", 
					  x = -0.1, 
                      y = 0.5, 
                      rot = 90, 
                      gp = gpar(fontsize = gbl.var$font.size, 
                                fontface = "bold"))
		} else {
			grid.text("-log(p-value)", 
                      x = -0.1, 
                      y = 0.5, 
                      rot = 90, 
                      gp = gpar(fontsize = gbl.var$font.size, 
                                fontface = "bold"))
		}
	
		#-------------------EQUIDISPOS BEGINS-------------
	
		gbl.var$equidis.pos <- NULL

		#DEBUG STATEMENT
		#cat("length(gbl.var$equidis.pos) ", length(gbl.var$equidis.pos), "\n")
		#cat("length(gbl.var$snp.num) ", length(gbl.var$snp.pos), "\n")

		gbl.var$snp.num <- length(gbl.var$sorted.snp.pos)

		if(config.var$DISP.LDMAP) {
			pos.increment <- (gbl.var$max.x - gbl.var$min.x) / gbl.var$snp.num
	  
			#Attempt to center the SNP labels 
			start.pos <- gbl.var$min.x
	  
			for(i in 0:(gbl.var$snp.num - 1)) {
				gbl.var$equidis.pos <- c(gbl.var$equidis.pos, pos.increment * i + start.pos)
			}
		} else {	
			pos.increment <- (gbl.var$max.x - gbl.var$min.x) / gbl.var$snp.num
	  
			#Attempt to center the SNP labels 
			start.pos <- gbl.var$min.x
	  		
			for(i in 0:(gbl.var$snp.num - 1)) {
				gbl.var$equidis.pos <- c(gbl.var$equidis.pos, pos.increment * i + start.pos)
			}
		}
		
		#-------------------CONNECTING LINES BEGINS-------------
		
		if(config.var$CONNECTING.LINES.ADJ == -1) {
			stop("CONNECTING.LINES.ADJ may not equal: -1")
		}

		tmp.correction.factor <- gbl.var$equidis.pos[1] * 10^(gbl.var$cur.exp - 9)	* (1 + config.var$CONNECTING.LINES.ADJ)	
                     
		#DEBUG STATEMENT
		#cat("tmp.correction.factor ", tmp.correction.factor, "\n")
		#cat("length(gbl.var$sorted.snp.pos) ", length(gbl.var$sorted.snp.pos), "\n")
		#cat("length(gbl.var$equidis.pos) ", length(gbl.var$equidis.pos), "\n")
		#cat("gbl.var$equidis.pos[1] ", gbl.var$equidis.pos[1], "\n")
		#cat("gbl.var$cur.exp ", gbl.var$cur.exp, "\n")
		#cat("config.var$CONNECTING.LINES.ADJ ", config.var$CONNECTING.LINES.ADJ, "\n")

		if(config.var$IMAGE.SIZE == 3.5) {
			y.finish.pos <- rep(-0.9, length(gbl.var$sorted.snp.pos)) * config.var$CONNECTING.LINES.FACTOR
		} else if(config.var$IMAGE.SIZE == 7) {
			y.finish.pos <- rep(-1.8, length(gbl.var$sorted.snp.pos)) * config.var$CONNECTING.LINES.FACTOR
		}

		#DEBUG STATEMENT
		#cat("gbl.var$axis.y[1] ", gbl.var$axis.y[1], "\n")
		#cat("y.finish.pos ", y.finish.pos, "\n")
		
		if(config.var$CONNECTING.LINES.VERT.ADJ != -1) {
			y.start.pos <- config.var$CONNECTING.LINES.VERT.ADJ
		} else if (config.var$IMAGE.SIZE == 3.5) {
			y.start.pos <- -0.5	
		} else if (config.var$IMAGE.SIZE == 7) {
			y.start.pos <- -0.7 
		} else {
			stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
		}
                                               
		if(config.var$DISP.LDMAP) {

			if(config.var$EVEN.SPACED) {
				if(config.var$IMAGE.SIZE == 3.5) {

					#cat("gbl.var$sorted.snp.pos ", gbl.var$sorted.snp.pos, "\n")
	
					connecting.lines <- segmentsGrob(x0 = (seq(1, length(gbl.var$sorted.snp.pos)) + config.var$CONNECTING.LINES.ADJ), 
													 x1 = seq(1, length(gbl.var$sorted.snp.pos)), 
													 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos, "char"), 
													 default.units = "native", 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else if(config.var$IMAGE.SIZE == 7) {
					connecting.lines <- segmentsGrob(x0 = (seq(1, length(gbl.var$sorted.snp.pos)) + config.var$CONNECTING.LINES.ADJ), 
													 x1 = seq(1, length(gbl.var$sorted.snp.pos)), 
													 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos,"char"), 
													 default.units = "native", 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else {
					stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
				}	
			} else {

				x.finish.pos <- gbl.var$equidis.pos + 
								seq(1, length(gbl.var$equidis.pos) * abs(tmp.correction.factor), abs(tmp.correction.factor)) * config.var$CONNECTING.LINES.FLEX + 
								gbl.var$total.dist * config.var$CONNECTING.LINES.ADJ

				if(config.var$IMAGE.SIZE == 3.5) {
					connecting.lines <- segmentsGrob(x0 = x.finish.pos, 
													 x1 = gbl.var$sorted.snp.pos, 
													 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos, "char"), 
													 default.units = "native", 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else if(config.var$IMAGE.SIZE == 7) {
					connecting.lines <- segmentsGrob(x0 = x.finish.pos, 
													 x1 = gbl.var$sorted.snp.pos, 
													 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos,"char"), 
													 default.units = "native", 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else {
					stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
				}	
			}	
		} else {
	
			if(config.var$EVEN.SPACED) {
				if(config.var$IMAGE.SIZE == 3.5) {
					connecting.lines <- segmentsGrob(x0 = (seq(1, length(gbl.var$sorted.snp.pos)) + config.var$CONNECTING.LINES.ADJ), 
													 x1 = seq(1, length(gbl.var$sorted.snp.pos)), 
	 											 	 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos, "char"), 
													 default.units = "native" , 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else if(config.var$IMAGE.SIZE == 7) {
					connecting.lines <- segmentsGrob(x0 = (seq(1, length(gbl.var$sorted.snp.pos)) + config.var$CONNECTING.LINES.ADJ), 
													 x1 = seq(1, length(gbl.var$sorted.snp.pos)), 
													 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos, "char"), 
													 default.units = "native" , 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else {
					stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
				}			
			} else {

				x.finish.pos <- gbl.var$equidis.pos + 
								seq(1, length(gbl.var$equidis.pos) * tmp.correction.factor, tmp.correction.factor) * config.var$CONNECTING.LINES.FLEX +
								gbl.var$total.dist * config.var$CONNECTING.LINES.ADJ

				if(config.var$IMAGE.SIZE == 3.5) {
					connecting.lines <- segmentsGrob(x0 = x.finish.pos, 
													 x1 = gbl.var$sorted.snp.pos, 
	 											 	 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos, "char"), 
													 default.units = "native" , 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else if(config.var$IMAGE.SIZE == 7) {
					connecting.lines <- segmentsGrob(x0 = x.finish.pos, 
													 x1 = gbl.var$sorted.snp.pos, 
													 y0 = unit(y.finish.pos, "char"), 
													 y1 = unit(y.start.pos, "char"), 
													 default.units = "native" , 
													 gp = gpar(lwd = gbl.var$line.width))	
				} else {
					stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
				}
			}
		}

		#DEBUG STATEMENT
		#cat("y.finish.pos ", y.finish.pos, "\n")
		#cat("connecting.lines ", is.null(connecting.lines), "\n")
		#cat(is.null(connecting.lines), "\n")
		#cat(x.finish.pos, "\n")
		#cat(gbl.var$sorted.snp.pos, "\n")

		if(config.var$DISP.CONNECTING.LINES) {
			grid.draw(connecting.lines)
		}
	
		#DEBUG STATEMENT
		cat("FINISH DRAW.PLOT.GRID.SETUP\n")

		return(gbl.var)
	}
	
	draw.legend <- function(config.var, gbl.var) {
	
		#DEBUG STATEMENT
		cat("START DRAW.LEGEND\n")

		#DEBUG STATEMENT
		#cat("gbl.var$font.size ", gbl.var$font.size, "\n")

		if(config.var$DISP.LDMAP) {
			legend.vp <- viewport(x = 0.5, y = -0.5, height = 1, width = 1, name = "legend.vp")
		} 
		else {
			legend.vp <- viewport(name = "legend.vp")
		}
		
		pushViewport(legend.vp)
	
		#-------------------LEGEND SAMPLES BEGINS-----------------------    
	
		legend.samples.list <- NULL
		legend.samples <- gTree(children=NULL, 
								 just=c("center", "bottom"), vp=legend.vp, 
								 name="legend.samples")
	
		if(samples <= 10) {
	
				for(i in 1:samples) {
	
					y.pos.sec.row.correction <- 0
					x.pos.sec.row.correction <- 0
	
					if(i > 5) {
						x.pos.sec.row.correction <- -1.5
						y.pos.sec.row.correction <- -.04
					}
	
					label.text <- substr(gbl.var$split.sample.labels[[1]][i], 1, 10)
	
					#Ensure that labels do not overlap by moving them
	       			if(config.var$DISP.LDMAP) {
						y.pos <- .75 + 0.05 * (i - 1)
						x.pos <- 0.85
	
					} else {
	                    y.pos <- -0.3125 + y.pos.sec.row.correction 
						x.pos <- -0.19 + 0.3 * (i - 1) + x.pos.sec.row.correction
					}
	
					label.sample.text <- textGrob(label.text, 
                                                   x = (x.pos + 0.035), 
                                                   y = y.pos, 
                                                   just=("left"), 
                                                   gp = gpar(fontsize = gbl.var$font.size))
	
					if(is.na(gbl.var$symbol.list[i])) {
						gbl.var$symbol.list[i] <- 24
					}	

					#FOR LD adjust is .3
					#FOR NO-LD adjust is 
	  				label.sample.symbol <- pointsGrob(x = x.pos, 
                                                       y = y.pos, 
                                                       pch = gbl.var$symbol.list[i],
										               gp = gpar(col = gbl.var$color.list[i], 
                                                                 cex = (gbl.var$cex.factor - 0.50), 
                                                                 lwd = gbl.var$line.width, 
                                                                 fill = gbl.var$fill.list[i]))
	
					legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.text
					legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.symbol
				}
	
			class(legend.samples.list) <- c("gList")
	
			legend.glist.samples <- legend.samples.list
	
			legend.samples <- gTree(children=legend.glist.samples, 
									 just=c("left", "bottom"), 
                                     vp=legend.vp, 
									 name="legend.samples")
	  
			grid.draw(legend.samples)
		}
	
		popViewport()
	
		#DEBUG STATEMENT
		cat("FINISH DRAW.LEGEND\n")
	}
	
	#draw SNP labels
	draw.plot.grid.snp.names <- function(config.var, gbl.var) {
	
		#DEBUG STATEMENT
		cat("START DRAW.PLOT.GRID.SNP.NAMES\n")

		tmp.correction.factor <- gbl.var$equidis.pos[1] * 10^(gbl.var$cur.exp - 9)	* (1 + config.var$CONNECTING.LINES.ADJ)

		#DEBUG STATEMENT
		#cat("config.var$DISP.LDMAP ", config.var$DISP.LDMAP, "\n")
		#cat("length(gbl.var$sorted.snp.names) ", length(gbl.var$sorted.snp.names), "\n")
		#cat("tmp.correction.factor ", tmp.correction.factor, "\n")
		#cat("gbl.var$axis.y[1] ", gbl.var$axis.y[1], "\n")
		
		if(!config.var$DISP.LDMAP) {

			if(config.var$EVEN.SPACED) {
				
				#DEBUG STATEMENT
				#cat("gbl.var$axis.y[1] ", gbl.var$axis.y[1], "\n")

				if(config.var$IMAGE.SIZE == 3.5) {
					snp.names.grob <- textGrob(gbl.var$sorted.snp.names, x = 1:length(gbl.var$sorted.snp.pos), y = unit(-8, "char"), 
										  rot=90, default.units = "native", 
										  gp=gpar(fontsize = gbl.var$font.size), just=c("left"),
										  name="snp.names")
				} else if(config.var$IMAGE.SIZE == 7) { 
					snp.names.grob <- textGrob(gbl.var$sorted.snp.names, x = 1:length(gbl.var$sorted.snp.pos), y = unit(-8, "char"), 
										  rot=90, default.units = "native", 
										  gp=gpar(fontsize = (gbl.var$font.size - 2)), just=c("left"),
										  name="snp.names")
				}
			} else { 		
				if(config.var$IMAGE.SIZE == 3.5) {
					snp.names.grob <- textGrob(gbl.var$sorted.snp.names, x = gbl.var$equidis.pos, y = unit(-10, "char"), 
										  rot=90, default.units = "native", 
										  gp=gpar(fontsize = (gbl.var$font.size - 1)), just=c("left"), name="snp.names")
				} else if (config.var$IMAGE.SIZE == 7) { 
					snp.names.grob <- textGrob(gbl.var$sorted.snp.names, x = gbl.var$equidis.pos, y = unit(-8, "char"), 
										  rot=90, default.units = "native", 
										  gp=gpar(fontsize = (gbl.var$font.size - 2)), just=c("left"), name="snp.names")
				}	
			}		

			if(config.var$DISP.SNP.NAMES) {
				#grid.draw(snp.names.grob)
			}
		}

		#DEBUG STATEMENT
		cat("FINISH DRAW.PLOT.GRID.SNP.NAMES\n")
	}
	
	#-------------------DRAW SINGLE SNPS------------------
	draw.plot.grid.ss <- function(config.var, gbl.var) {
	
		#DEBUG STATEMENT
		cat("START DRAW.PLOT.GRID.SS\n")

		# Grab from configuration is available 
		if(is.null(config.var$SYMBOL.FACTOR)) {
			tmp.cex.factor.symbol = gbl.var$cex.factor.symbol
		} else {
			tmp.cex.factor.symbol = config.var$SYMBOL.FACTOR
		}
	
		if(config.var$DISP.SNP) {
		  	if(config.var$DISP.TYPE == "allele") {
		  		grid.text(gbl.var$snp.data$MAJOR.ALLELE, gbl.var$sorted.snp.pos, gbl.var$snp.data$SS.PVAL, 
						  default.units = "native", gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], cex = tmp.cex.factor.symbol))
		  	}
		  	else if(config.var$DISP.TYPE == "base") {
		  		for (i in 1:nrow(gbl.var$snp.data)) {
		  
		  			if(gbl.var$snp.data$SS.PVAL[i] > config.var$PVAL.THRESHOLD) {
		  				if(gbl.var$snp.data$MAJOR.ALLELE[i] == 1) {
		  					gbl.var$snp.data$MAJOR.ALLELE[i] <- gbl.var$snp.bases[[i]][1]
		  					grid.text(gbl.var$snp.data$MAJOR.ALLELE[i], 
									  gbl.var$snp.data$LOC[i], 
									  gbl.var$snp.data$SS.PVAL[i], 
									  default.units = "native", 
									  gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample],
									  			fontface = "italic", 
												cex = tmp.cex.factor.symbol))
		  				}
		  				else if(gbl.var$snp.data$MAJOR.ALLELE[i] == 2) {
		  					gbl.var$snp.data$MAJOR.ALLELE[i] <- gbl.var$snp.bases[[i]][2]
		  					grid.text(gbl.var$snp.data$MAJOR.ALLELE[i], 
									  gbl.var$snp.data$LOC[i], 
									  gbl.var$snp.data$SS.PVAL[i], 
									  default.units = "native", 
									  gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
									  cex = tmp.cex.factor.symbol))
		  				}
		  			}
		  		}
	  		}	
		  	else if(config.var$DISP.TYPE == "symbol") {
		  
		  		for (i in 1:nrow(gbl.var$snp.data)) {
		  			if ((gbl.var$snp.data$SS.PVAL[i] > config.var$PVAL.THRESHOLD) & (gbl.var$snp.data$SS.PVAL[i] != 0)) {
		  
		  				if (is.na(gbl.var$symbol.list[gbl.var$cur.sample])) {
		  					if (gbl.var$snp.data$ASSOC[i] == "+") {
		  						grid.points(gbl.var$snp.data$LOC[i], 
											gbl.var$snp.data$SS.PVAL[i], 
											pch = 24, 
											gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
													  cex = tmp.cex.factor.symbol, 
													  fill = gbl.var$fill.list[gbl.var$cur.sample],
													  lwd = gbl.var$line.width))
		  					}
		  					else if (gbl.var$snp.data$ASSOC[i] == "-") {
		  						grid.points(gbl.var$snp.data$LOC[i], 
											gbl.var$snp.data$SS.PVAL[i], 
											pch = 25, 
											gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
													  cex = tmp.cex.factor.symbol,
													  fill = gbl.var$fill.list[gbl.var$cur.sample], 
													  lwd = gbl.var$line.width))
		  					}
		  				}
		  				else { 
		  					grid.points(gbl.var$snp.data$LOC[i], 
										gbl.var$snp.data$SS.PVAL[i], 
										pch = gbl.var$symbol.list[gbl.var$cur.sample],
										gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
												  cex = tmp.cex.factor.symbol, 
												  fill = gbl.var$fill.list[gbl.var$cur.sample], 
												  lwd = gbl.var$line.width))
		  				}
		  			}
		  		}	
	  		}
	  		else {
	  			stop("Invalid display type: ", config.var$DISP.TYPE, "\n")
	  		}
		}
	
		#DEBUG STATEMENT
		cat("FINISH DRAW.PLOT.GRID.SS\n")
	}
	
	#-------------------HAPLOTYPES BEGINS------------------
	draw.plot.grid.hap <- function(config.var, gbl.var) { 
	
		#DEBUG STATEMENT
		cat("START DRAW.PLOT.GRID.HAP\n")

		tmp.colnames.loc <- NULL
		
		# Grab from configuration is available 
		if(is.null(config.var$SYMBOL.FACTOR)) {
			tmp.cex.factor.symbol = gbl.var$cex.factor.symbol
		} else {
			tmp.cex.factor.symbol = config.var$SYMBOL.FACTOR
		}
	
		if(config.var$DISP.HAP) {
	
			#Indicates the position of the first SNP
			hap.offset <- 4
	
			#if(length(hap.offset) == 0) {
			#	stop("SNP name mismatch between the SNP and haplotype data files in the sample: ", gbl.var$cur.sample)
			#}
	
			if(config.var$USE.GBL.PVAL) {
				for(i in 1:nrow(gbl.var$hap.data)) {
	
					#Look for the haplotype in each row if it has a significant global p-value
					if(!is.na(gbl.var$hap.data$G.PVAL[i]) & (gbl.var$hap.data$G.PVAL[i] > config.var$PVAL.THRESHOLD)) {
	
						#cur.haplo <- which(!is.na(gbl.var$hap.data[i,hap.offset:ncol(gbl.var$hap.data)])) + (hap.offset - 1)
						cur.haplo <- which(!is.na(gbl.var$hap.data[i,hap.offset:ncol(gbl.var$hap.data)])) + (hap.offset - 1)
						cur.haplo.colnames <- colnames(gbl.var$hap.data)[cur.haplo]
	
						for(j in 1:length(cur.haplo.colnames)) {
							#DEBUG STATEMENT
							#cat("cur.haplo.colnames[j], ", cur.haplo.colnames[j], "\n")
							#cat("grep ", grep(paste(cur.haplo.colnames[j],  "$", sep=""), gbl.var$snp.data$SNP.NAME, ignore.case = TRUE, value = TRUE), "\n")

							tmp.colnames.loc <- c(tmp.colnames.loc, grep(paste(cur.haplo.colnames[j],  "$", sep=""), gbl.var$snp.data$SNP.NAME, ignore.case = TRUE))
						}

						#DEBUG STATEMENT
						#cat("TCL: ", tmp.colnames.loc, "\n")
						#cat("CH: ", cur.haplo, "\n")
	
						if(length(gbl.var$snp.data$LOC[tmp.colnames.loc]) != length(cur.haplo)) {
							stop("Please make sure each column in the haplotype data file correctly corresponds to a SNP in the SNP data file.\n")
						}

						#grid.text(gbl.var$hap.data[i, cur.haplo], gbl.var$snp.data$LOC[cur.haplo - (hap.offset - 1)], rep(gbl.var$hap.data$G.PVAL[i], length(cur.haplo)))
						grid.text(gbl.var$hap.data[i, cur.haplo], 
								  gbl.var$snp.data$LOC[tmp.colnames.loc], 
								  rep(gbl.var$hap.data$G.PVAL[i], length(cur.haplo)))	
						
						if(is.na(gbl.var$symbol.list[gbl.var$cur.sample])) {	
							#cat("SIZE ", length(gbl.var$snp.data$LOC[tmp.colnames.loc]), "\n")
							grid.points(gbl.var$snp.data$LOC[tmp.colnames.loc], 
										rep(gbl.var$hap.data$G.PVAL[i], length(cur.haplo)), 
										pch = 21, 
										gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
												  cex = tmp.cex.factor.symbol, 
												  lwd = gbl.var$line.width))
						} else {
							grid.points(gbl.var$snp.data$LOC[tmp.colnames.loc], 
										rep(gbl.var$hap.data$G.PVAL[i], length(cur.haplo)), 
										pch = gbl.var$symbol.list[gbl.var$cur.sample], 
										gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
												  cex = tmp.cex.factor.symbol, 
												  lwd = gbl.var$line.width, 
												  fill = gbl.var$fill.list[gbl.var$cur.sample]))
						}
	
						grid.lines(gbl.var$snp.data$LOC[tmp.colnames.loc], 
								   rep(gbl.var$hap.data$G.PVAL[i], length(cur.haplo)), 
								   default.units = "native", 
								   gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
											 lwd = gbl.var$line.width))
	
						tmp.pval <- gbl.var$hap.data$G.PVAL[i]
						tmp.colnames.loc <- NULL
					}
				}
			}
	
			# Draw haplotyes, if positive association display solid line, if negative display dashed line. 
			if(!config.var$USE.GBL.PVAL) {
				for(i in 1:nrow(gbl.var$hap.data)) {
		
					#DEBUG STATEMENT 
					#cat("gbl.var$hap.data$I.PVAL[i] ", gbl.var$hap.data$I.PVAL[i], "\n")
					#cat("config.var$PVAL.THRESHOLD ", config.var$PVAL.THRESHOLD, "\n")					

					#The ln() of the p-values is being checked
					if(!is.na(gbl.var$hap.data$I.PVAL[i]) & (gbl.var$hap.data$I.PVAL[i] > config.var$PVAL.THRESHOLD)) {
	
						cur.haplo <- which(!is.na(gbl.var$hap.data[i,hap.offset:ncol(gbl.var$hap.data)])) + (hap.offset - 1)
						cur.haplo.colnames <- colnames(gbl.var$hap.data)[cur.haplo]
	
						for(j in 1:length(cur.haplo.colnames)) {
							#DEBUG STATEMENT
							#cat("cur.haplo.colnames[j], ", cur.haplo.colnames[j], "\n")
							#cat("grep ", grep(paste(cur.haplo.colnames[j],  "$", sep=""), gbl.var$snp.data$SNP.NAME, ignore.case = TRUE, value = TRUE), "\n")

							tmp.colnames.loc <- c(tmp.colnames.loc, grep(paste(cur.haplo.colnames[j],  "$", sep=""), gbl.var$snp.data$SNP.NAME, ignore.case = TRUE))
						}

						#DEBUG STATEMENT
						#cat("TCL: ", tmp.colnames.loc, "\n")
						#cat("CH: ", cur.haplo, "\n")

						if(length(gbl.var$snp.data$LOC[tmp.colnames.loc]) != length(cur.haplo)) {
							stop("Please make sure each column in the haplotype data file correctly corresponds to a SNP in the SNP data file.\n")
						}
	
						if(config.var$DISP.TYPE == "allele") {
							grid.text(gbl.var$hap.data[i, cur.haplo], 
									  gbl.var$snp.data$LOC[tmp.colnames.loc], 
									  rep(gbl.var$hap.data$I.PVAL[i], length(cur.haplo)), 
									  default.units = "native", 
									  gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
												cex = gbl.var$cex.factor))
						} else if(config.var$DISP.TYPE == "base") {
							for (j in tmp.colnames.loc) {
							
								if(gbl.var$hap.data[i, j] == 1) {
									#gbl.var$hap.data[i, j] <- gbl.var$snp.bases[[j - (hap.offset - 1)]][1]
									gbl.var$hap.data[i, j] <- gbl.var$snp.bases[[j]][1]
									
									grid.text(gbl.var$hap.data[i, j], 
											  gbl.var$snp.data$LOC[j - (hap.offset - 1)], 
											  gbl.var$hap.data$I.PVAL[i], 
											  default.units = "native", 
											  gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
														fontface = "italic", 
														cex = gbl.var$cex.factor))
								} else if(gbl.var$hap.data[i, j] == 2) {
									gbl.var$hap.data[i, j] <- gbl.var$snp.bases[[j]][2]
									grid.text(gbl.var$hap.data[i, j], 
											  gbl.var$snp.data$LOC[j - (hap.offset - 1)], 
											  gbl.var$hap.data$I.PVAL[i], 
											  default.units = "native", 
											  gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
														cex = gbl.var$cex.factor))
								} else {
									stop("Only DISP.TYPE 'allele' and 'symbol' can display more than 2 alleles.\n")
								}
							}
						} else if(config.var$DISP.TYPE == "symbol") {
								for (j in 1:length(cur.haplo)) {

									if (is.na(gbl.var$symbol.list[gbl.var$cur.sample])) {
										gbl.var$symbol.list[gbl.var$cur.sample] <- 21
									}

									if(gbl.var$hap.data[i, cur.haplo[j]] == 1) {
										grid.points(gbl.var$snp.data$LOC[tmp.colnames.loc[j]], 
													gbl.var$hap.data$I.PVAL[i], 
													pch = gbl.var$symbol.list[gbl.var$cur.sample], 
													gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
															  cex = tmp.cex.factor.symbol, 
															  lwd = gbl.var$line.width, 
															  fill = gbl.var$color.list[gbl.var$cur.sample]))
									} else if(gbl.var$hap.data[i, cur.haplo[j]] == 2) {
										grid.points(gbl.var$snp.data$LOC[tmp.colnames.loc[j]], 
													gbl.var$hap.data$I.PVAL[i], 
													pch = gbl.var$symbol.list[gbl.var$cur.sample], 
													gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
															  cex = tmp.cex.factor.symbol, 
															  lwd = gbl.var$line.width))
									}
									else {
										stop("Only DISP.TYPE 'allele' and 'symbol' can display more than 2 alleles.\n")
									}
								}
						}	
						
						#Draw lines between haplotype points based on association 
						#If NA is given as association, then draw a solid line
						if(gbl.var$hap.data$ASSOC[i] == "+") {
							grid.lines(gbl.var$snp.data$LOC[tmp.colnames.loc], 
									   rep(gbl.var$hap.data$I.PVAL[i], length(cur.haplo)), 
									   default.units = "native", 
									   gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
												 cex = gbl.var$cex.factor, 
												 lwd = (gbl.var$line.width)))
						} else if(gbl.var$hap.data$ASSOC[i] == "-") {
							grid.lines(gbl.var$snp.data$LOC[tmp.colnames.loc], 
									   rep(gbl.var$hap.data$I.PVAL[i], length(cur.haplo)), 
									   default.units = "native", 
									   gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
                                                 cex = gbl.var$cex.factor, 
                                                 lwd = (gbl.var$line.width), 
                                                 lty = "dashed"))
						} else if(is.na(gbl.var$hap.data$ASSOC[i])) {
							grid.lines(gbl.var$snp.data$LOC[tmp.colnames.loc], 
									   rep(gbl.var$hap.data$I.PVAL[i], length(cur.haplo)), 
                                       default.units = "native", 
                                       gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
                                                 cex = gbl.var$cex.factor, 
                                                 lwd = gbl.var$line.width))
						} else {
							stop("Association could not be found for 1 or more SNPs. A generic symbol SYMBOLS must be specified\n")
						}	

						if(length(tmp.colnames.loc) >= length(cur.haplo)) {
							tmp.colnames.loc <- NULL
						}	
					}
				}
			}
		}
	
		#DEBUG STATEMENT
		cat("FINISH DRAW.PLOT.GRID.HAP\n")
	}
	
	draw.plot.ld.plot <- function(config.var, gbl.var) { 
	
		#DEBUG STATEMENT
		cat("START DRAW.PLOT.LD.PLOT\n")
	
		sample.ld <- LD(gbl.var$geno.data)

		if(config.var$LD.TYPE == "rsquare") {
			ld.matrix <- sample.ld[["r"]]^2
		} else if(config.var$LD.TYPE == "dprime") {
			ld.matrix <- abs(sample.ld[["D'"]])
		} else {
			stop("Invalid LD type: ", config.var$LD.TYPE, "\n")
		}
	
		if(config.var$USE.COLORS) {
			if(config.var$LD.COLOR.SCHEME == "heat") {
				ld.colors <- heat.colors(gbl.var$palette.size)
			} else if(config.var$LD.COLOR.SCHEME == "cm") {
					ld.colors <- cm.colors(gbl.var$palette.size)
			} else if(config.var$LD.COLOR.SCHEME == "custom") {
				if(!is.null(config.var$PALETTE.FILE)) {
						tmp.colors <- read.table(config.var$PALETTE.FILE, as.is=TRUE, header=FALSE, blank.lines.skip = TRUE)
						gbl.var$palette.size <- length(tmp.colors$V1)
						custom.colors <- paste("#", tmp.colors$V1, sep="")
						ld.colors <- rev(custom.colors)
				}
				else {
					stop("PALETTE.FILE must be specified\n")
				}
			}
			else if(config.var$LD.COLOR.SCHEME == "topo") {
				ld.colors <- topo.colors(gbl.var$palette.size)
			}
			else if(config.var$LD.COLOR.SCHEME == "gray") {
				ld.colors <- gray.colors(gbl.var$palette.size)
			}
			else {
				stop("Invalid color scheme: ", config.var$LD.COLOR.SCHEME, "\n")
			}
		}
		else {
	        ld.colors <- gray.colors(gbl.var$palette.size)
		}

		color.bar.colors <- c(rep(NA, gbl.var$palette.size), ld.colors[gbl.var$palette.size:1])	

		#DEBUG STATEMENT 
		#cat("color.bar.colors ", color.bar.colors, "\n")

		color.bar <- rectGrob(x=rep(seq(1:gbl.var$palette.size)/gbl.var$palette.size, 2), 
							  y=rep(seq(1:2)/2, each=gbl.var$palette.size), 
                              width=(1/gbl.var$palette.size), 
                              height=0.125, 
                              just=c("right", "top"), 
                              gp=gpar(col=NULL, 
                                      fill=color.bar.colors, 
                                      cex = (gbl.var$cex.factor + 0.2), 
                                      lty="blank"), 
                              name = "color.bar")
	
	  	color.bar.labels <- textGrob(paste(0.2*0:5), x=0.2*0:5, y=0.8, gp=gpar(fontsize = gbl.var$font.size), name="color.bar.labels")
	
		ld.key.vp <- viewport(x=0.2, y=-0.03, width = 0.25, height = 0.25)
		ld.key <- gTree(children=gList(color.bar, color.bar.labels), name = "ld.key", vp=ld.key.vp)	
	
		if(config.var$DISP.PHYS.DIST) {
			if(gbl.var$total.dist > 1000) {
				map.label.text.distance <- paste("Physical Distance: ", round(gbl.var$total.dist/1000, 1), " kb", sep="")
			} else {
				map.label.text.distance <- paste("Physical Distance: ", gbl.var$total.dist, " bases", sep="")
			}
		} else {
			map.label.text.distance <- " "
		}
		
		if (config.var$LD.TYPE == "rsquare") {
			map.label.text.ldtype <- "LD Map Type: r-square"
		} else if (config.var$LD.TYPE == "dprime") {
			map.label.text.ldtype <- "LD Map Type: D'"
		} else {
			stop("Invalid LD metric : ", config.var$LD.TYPE, "\n")
		}
	
		map.label.distance <- textGrob(map.label.text.distance, 
									   x = 0.2, 
									   y = 0.1425, 
									   just=("center"), 
									   gp = gpar(fontsize = gbl.var$font.size), 
											     name="map.label.distance")
		map.label.ldtype <- textGrob(map.label.text.ldtype, 
									 x = 0.2, 
									 y = 0.11150, 
									 just=("center"), 
									 gp = gpar(fontsize = gbl.var$font.size), 
									 name="map.label.ldtype")
	
		matrix.rows <- dim(ld.matrix)[1]
		matrix.cols <- dim(ld.matrix)[2] 
	
		color.intervals <- 0:length(ld.colors)/length(ld.colors)
		tmp.color.cut <- as.character(cut(1-ld.matrix, color.intervals, labels=as.character(ld.colors)))
	
		tmp.vector <- NULL
	
		for(i in seq(0, matrix.rows - 1)) {
			for(j in seq(matrix.rows, 1, by =-1)) {
				tmp.vector <- c(tmp.vector, (matrix.rows*j-i))			
			}
		}
		
		color.cut <- rep(NULL, matrix.rows^2)
	
		for(i in 1:matrix.rows^2) {
			color.cut[i] <- tmp.color.cut[tmp.vector[i]]
		}
	
		x.x <- (1:matrix.cols)/matrix.cols
		y.y <- (1:matrix.rows)/matrix.rows
		right <- rep(x.x, nrow(gbl.var$geno.data))
		top <- rep(y.y, each=matrix.cols)
	
		image.rect <- rectGrob(x=right, 
							   y=top, 
                               width=1/matrix.cols, 
                               height=1/matrix.rows, 
                               just=c("right", "top"), 
                               gp=gpar(col=NULL, 
                                       fill=color.cut, 
                                       lty="blank"), 
                                       name="image.rect")
	
		#seq.x <- c(1*(1/gbl.var$snp.num), 1)
		#seq.y <- c(0, (1 - (1/gbl.var$snp.num)))
	
		tmp.snp.num <- gbl.var$snp.num
	
		seq.x.names <- seq((1/tmp.snp.num), 1, by=(1/tmp.snp.num))
		seq.y.names <- seq(0, (1 - (1/tmp.snp.num)), by=(1/tmp.snp.num))
	
		#The condition for the tmp.snp.num is met in the original construction of the sequence
		if(tmp.snp.num > 50) {
			tmp.x.factor <- seq.x.names[1]		
	
			seq.x.names <- seq.x.names + tmp.x.factor/2
			seq.y.names <- seq.y.names - tmp.x.factor/2
		} else if(tmp.snp.num < 50) {
			tmp.x.factor <- (seq.x.names[1] - 0.02)/2 + 0.02
			tmp.y.factor <- (seq.x.names[1] - 0.02)/2
	
			seq.x.names <- seq(tmp.x.factor, (1 - tmp.y.factor), by=(1/tmp.snp.num))
			seq.y.names <- seq(tmp.y.factor, (1 - tmp.y.factor), by=(1/tmp.snp.num))
		} 
	
		snp.names.pos.x <- (1:gbl.var$snp.num)/gbl.var$snp.num
		snp.names.pos.y <- ((1:gbl.var$snp.num)/gbl.var$snp.num - 1/gbl.var$snp.num)

		if(is.null(config.var$FONT.FACTOR)) {
			tmp.font.factor = (gbl.var$cex.factor - 0.5)
		} else {
			tmp.font.factor = config.var$FONT.FACTOR
		}

		if (config.var$DISP.SNP.NAMES) {						
			snp.names.ld.sec <- textGrob(rev(gbl.var$sorted.snp.names), 
										 x = seq.x.names[1:gbl.var$snp.num], 
										 y = seq.y.names[1:gbl.var$snp.num], 
										 rot=-45, 
										 gp=gpar(cex = tmp.font.factor), 
										 just =  c("left"), 
										 name="snp.names.ld.sec")
		}

	  	popViewport()
	
		#67.5
		ld.map.vp <- viewport(y = unit(0.415, "npc"), 
							  x = unit(0.471875, "npc"), 
							  width = unit(0.7, "snpc"), 
							  height = unit(0.7, "snpc"), 
							  angle = 67.5, 
							  just = c("center", "center"), 
							  name = "ld.map.vp")
		if (config.var$DISP.SNP.NAMES) {
			ld.map <- gTree(children=gList(image.rect, snp.names.ld.sec), just=c("center", "bottom"), vp=ld.map.vp, name="ld.map")	
		} else {
			ld.map <- gTree(children=gList(image.rect), just=c("center"), vp=ld.map.vp, name="ld.map")			
		}
		
		gene.map.vp <- viewport(name = "gene.map.vp") 
		gene.map <- gTree(children=gList(map.label.ldtype, map.label.distance), vp = gene.map.vp, name="gene.map")
		
		pushViewport(ld.map.vp)
	
		grid.draw(ld.map) 
	
		popViewport() 
	
		pushViewport(gene.map.vp) 
	
		grid.draw(gene.map) 

		if(config.var$DISP.COLOR.BAR) {
			#DEBUG STATEMENT
			#cat("is.null(ld.key) ", is.null(ld.key), "\n")

			grid.draw(ld.key)
		}
		
		popViewport()
	
		#DEBUG STATEMENT
		cat("FINISH DRAW.PLOT.LD.PLOT\n")
	}
	
	set.image.parameters <- function(config.var, gbl.var) {
	
		#DEBUG STATEMENT
		cat("START SET.IMAGE.PARAMETERS\n")
	
		sec.cex.factor <- 1
	
		#DEBUG STATEMENT
		#cat("IMAGE.SIZE ", config.var$IMAGE.SIZE, "\n")

		if(config.var$IMAGE.SIZE == 3.5) {
			font.size <- 5
			line.width <- 0.5
			sec.cex.factor <- 1.5
			cex.factor.symbol <- 0.25
		} else if(config.var$IMAGE.SIZE == 7) {
			font.size <- 12
			line.width <- 1
			sec.cex.factor <- 2
			cex.factor.symbol <- 0.5
		} else {
			stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
		}
	
		if(gbl.var$snp.num > 0 & gbl.var$snp.num <= 20 ) {
			cex.factor <- 0.54 * sec.cex.factor
		} else if(gbl.var$snp.num > 20 & gbl.var$snp.num <= 60 ) {
			cex.factor <- 0.54 * sec.cex.factor
		} else if(gbl.var$snp.num > 60 & gbl.var$snp.num <= 90) {
			cex.factor <- 0.4 * sec.cex.factor
		} else if(gbl.var$snp.num > 90 & gbl.var$snp.num <= 120) {
			cex.factor <- 0.2 * sec.cex.factor
		} else if(gbl.var$snp.num > 150) {
			cex.factor <- 0.1 * sec.cex.factor
			line.width <- 0.5
		}

		#DEBUG STATEMENT
		#cat("gbl.var$snp.num ", gbl.var$snp.num, "\n")
		#cat("font.size ", font.size, "\n")
		#cat("line.width ", line.width, "\n")
		#cat("cex.factor.symbol ", cex.factor.symbol, "\n")
		#cat("cex.factor ", cex.factor, "\n")
		#cat("sec.cex.factor ", sec.cex.factor, "\n")
	
		#DEBUG STATEMENT
		cat("FINISH SET.IMAGE.PARAMETERS\n")

		return(list(font.size = font.size, line.width = line.width, cex.factor.symbol = cex.factor.symbol, cex.factor = cex.factor))
	}
	
	create.color.list <- function(config.var) {
	
		#DEBUG STATEMENT
		cat("START CREATE.COLOR.LIST\n")
	
		if(config.var$USE.COLORS) {
			tmp.color.list <-  c("red", "blue", "green", "black", "orange")
		}
		else {
	        tmp.color.list <-  c("grey0", "grey20", "grey40", "grey60", "grey80")
		}
	
		if(!is.null(config.var$COLOR.LIST)) {
			split.tmp.color.list <- strsplit(config.var$COLOR.LIST, ",")
	
			if(gbl.var$snp.samples != length(split.tmp.color.list[[1]])) { 
				stop("The color list must have ", samples, " colors separated by commas without spaces.\n")
			}
		} else {
			if(gbl.var$snp.samples <= 5) {
				split.tmp.color.list <- list(head(tmp.color.list, gbl.var$snp.samples))
			} else {
				stop("snp.plotter includes a default set of colors for 5 samples. For more samples, please specify COLOR.LIST\n")
	  		}
		}
	
		#DEBUG STATEMENT
		cat("FINISH CREATE.COLOR.LIST\n")
	
		return(split.tmp.color.list)
	}
	
	create.symbol.list <- function(config.var, split.color.list, gbl.var) {
	
		#DEBUG STATEMENT
		cat("START CREATE.SYMBOL.LIST\n")
	
		tmp.symbol.list <- c("circle-fill", "square-fill", "diamond-fill", "triangle-fill", "circle")
	
		if(!is.na(config.var$SYMBOLS)) {
			split.tmp.symbol.list <- strsplit(config.var$SYMBOLS, ",")
	
	        if(gbl.var$snp.samples != length(split.tmp.symbol.list[[1]])) {
				stop("The symbol list must have ", gbl.var$snp.samples, " symbol(s) separated by commas without spaces.\n")
			}
		} else {
			if(gbl.var$snp.samples <= 5) {
	            split.tmp.symbol.list <- list(head(tmp.symbol.list, samples))
			} else {
	   			stop("snp.plotter includes a default set of colors for 5 samples. For more samples, please specify SYMBOLS\n")
			}
		}
	
		split.symbol.list <- NULL
	
		for(i in 1:gbl.var$snp.samples) {
			if(any(grep("circle", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
				split.symbol.list <- c(split.symbol.list, 21)
			} else if(any(grep("square", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
	            split.symbol.list <- c(split.symbol.list, 22)
			} else if(any(grep("diamond", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
	            split.symbol.list <- c(split.symbol.list, 23)
			} else if(any(grep("triangle", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
	            split.symbol.list <- c(split.symbol.list, 24)
			} else if(any(grep("NA", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
	            split.symbol.list <- c(split.symbol.list, NA)
			} else {
				stop("Unknown symbol: ", split.tmp.symbol.list[[1]][i], "\n")
			}
		}
	
		split.fill.list <- list(rep("white", gbl.var$snp.samples))
	
		fill.pos <- grep("fill", split.tmp.symbol.list[[1]], ignore.case = TRUE)
		
		#DEBUG STATEMENT
		#cat("fill.pos ", fill.pos, "\n")

	   	split.fill.list[[1]][fill.pos] <- split.color.list[[1]][fill.pos]
	
		#DEBUG STATEMENT
		cat("FINISH CREATE.SYMBOL.LIST\n")
	
		return(list(split.symbol.list = split.symbol.list, split.fill.list = split.fill.list))
	}

	#-------------------GLOBAL VARIABLES BEGINS------------------

	#create three variables (global variables), which are lists of all other variables 
	#used throughout the program to be passed to functions with values changed as
	#necessary. 
	gbl.var <- NULL

	#DATA VARIABLES
	
	snp.data <- NULL
	hap.data <- NULL
	geno.data <- NULL

	snp.headers <- NULL
	snp.bases <- NULL
	snp.num <- 0
	snp.names <- NULL
	
	sorted.snp.names <- NULL
	sorted.snp.pos <- NULL
	
	cur.sample <- NULL
	snp.samples <- NULL
	hap.samples <- NULL

	snp.hash.names.pos <- new.env(hash=T)
	snp.hash.pos.names <- new.env(hash=T)
		
	#GRAPH BOUNDARY VARIABLES

	pval.flag <- TRUE 

	cur.exp <- 0

	total.dist <- NULL
	
	min.x <- NULL
	min.y <- NULL
	 
	max.x <- NULL
	max.y <- NULL
	
	r.x <- NULL
	r.y <- NULL
	
	axis.x <- NULL
	axis.y <- NULL

	equidis.pos <- NULL

	#FORMATTING VARIABLES
	
	split.sample.labels <- NULL
	split.color.list <- NULL
	color.list <- NULL
	symbol.list <- NULL
	fill.list <- NULL	

	palette.size <- 20
	cex.factor <- 0.5
	cex.factor.symbol <- 0.25
	font.size <- NULL
	line.width <- 0.5

	gbl.var <- list(snp.data = snp.data, 
					 hap.data = hap.data, 
					 geno.data = geno.data, 
					 snp.bases = snp.bases,
					 snp.num = snp.num,
					 snp.names = snp.names,
					 sorted.snp.names = sorted.snp.names,
					 sorted.snp.pos = sorted.snp.pos,
					 cur.sample = cur.sample, 
					 snp.samples = snp.samples, 
					 hap.samples = hap.samples,
 					 snp.hash.names.pos = snp.hash.names.pos, 
					 snp.hash.pos.names = snp.hash.pos.names,
				 	 pval.flag = pval.flag, 
					 cur.exp = cur.exp,
					 total.dist = total.dist,
					 min.x = min.x, 
					 min.y = min.y, 
					 max.x = max.x, 
					 max.y = max.y, 
					 r.x = r.x, 
				 	 r.y = r.y,
					 axis.x = axis.x, 
					 axis.y = axis.y, 
					 equidis.pos = equidis.pos, 
					 split.sample.labels = split.sample.labels,
					 split.color.list = split.color.list, 
				 	 color.list = color.list, 
					 symbol.list = symbol.list,
					 fill.list = fill.list,
					 palette.size = palette.size, 
					 cex.factor = cex.factor,
					 cex.factor.symbol = cex.factor.symbol, 
					 font.size = font.size, 
					 line.width = line.width
					)

	#-------------------GLOBAL VARIABLES ENDS------------------	

	#-------------------CONFIGURATION VARIABLES BEGINS---------
	
	config.var <- list(SNP.FILE = SNP.FILE, 
			 HAP.FILE = HAP.FILE, 
			 PALETTE.FILE = PALETTE.FILE, 
			 LAB.Y = LAB.Y, 
			 SYMBOLS = SYMBOLS,
			 EVEN.SPACED = EVEN.SPACED, 
			 PVAL.THRESHOLD = PVAL.THRESHOLD, 
			 USE.GBL.PVAL = USE.GBL.PVAL, 
			 DISP.TYPE = DISP.TYPE, 
			 DISP.LDMAP = DISP.LDMAP, 
			 DISP.HAP = DISP.HAP, 
			 DISP.SNP = DISP.SNP, 
			 DISP.MARKER.LINES = DISP.MARKER.LINES, 
			 DISP.CONNECTING.LINES = DISP.CONNECTING.LINES, 
			 DISP.SNP.NAMES = DISP.SNP.NAMES, 			 
			 LD.COLOR.SCHEME = LD.COLOR.SCHEME, 
			 COLOR.LIST = COLOR.LIST, 
			 USE.COLORS = USE.COLORS,
			 LD.TYPE = LD.TYPE, 
			 DISP.COLOR.BAR = DISP.COLOR.BAR,
			 DISP.PHYS.DIST = DISP.PHYS.DIST,
			 GENOTYPE.FILE = GENOTYPE.FILE, 
		 	 IMAGE.TITLE = IMAGE.TITLE, 
			 DISP.LEGEND = DISP.LEGEND, 
			 SAMPLE.LABELS = SAMPLE.LABELS, 
			 IMAGE.TYPE = IMAGE.TYPE, 
			 IMAGE.SIZE = IMAGE.SIZE,
			 DISP.MULT.LAB.X = DISP.MULT.LAB.X, 
			 IMAGE.NAME = IMAGE.NAME,
			 CONNECTING.LINES.FACTOR = CONNECTING.LINES.FACTOR,
			 CONNECTING.LINES.ADJ = CONNECTING.LINES.ADJ,
			 CONNECTING.LINES.VERT.ADJ = CONNECTING.LINES.VERT.ADJ,
			 CONNECTING.LINES.FLEX = CONNECTING.LINES.FLEX,
			 FONT.FACTOR = FONT.FACTOR,
			 SYMBOL.FACTOR = SYMBOL.FACTOR)

	#-------------------CONFIGURATION VARIABLES BEGINS---------
	
	if(!is.null(config.file)) {
		config.var <- read.config(config.file, config.var)
	}
	
	#START IMAGE CAPTURE
	if(config.var$IMAGE.TYPE == "pdf") {

		tmp.image.name <- paste(config.var$IMAGE.NAME, ".pdf", sep="")
		pdf(encoding = "ISOLatin1.enc", 
                   file = tmp.image.name, 
                   onefile=FALSE, 
                   width=as.numeric(config.var$IMAGE.SIZE),
                   height=as.numeric(config.var$IMAGE.SIZE), 
                   paper="special")
	} else if(config.var$IMAGE.TYPE == "eps") {
	
		tmp.image.name <- paste(config.var$IMAGE.NAME, ".eps", sep="")
		postscript(encoding = "ISOLatin1.enc", 
                   file = tmp.image.name, 
                   horizontal=FALSE, 
                   onefile=FALSE, 
                   width=as.numeric(config.var$IMAGE.SIZE), 
                   height=as.numeric(config.var$IMAGE.SIZE), 
                   paper="special", 
                   pagecentre=TRUE)
	} else {
		stop("Invalid image type: ", config.var$IMAGE.TYPE, " and image size: ", config.var$IMAGE.SIZE, " combination.\n")
	}
	#END IMAGE CAPTURE
	
	split.snp.file <- strsplit(config.var$SNP.FILE, ",")
	
	if(config.var$DISP.HAP) {
		split.hap.file <- strsplit(config.var$HAP.FILE, ",")
	}
	
	snp.file.length <- length(split.snp.file[[1]])
	gbl.var$snp.samples <- snp.file.length
    	gbl.var$hap.samples <- 0
	samples <- snp.file.length

	if(config.var$DISP.HAP) {
		split.hap.file <- strsplit(config.var$HAP.FILE, ",")
		hap.file.length <- length(split.hap.file[[1]])
	  
		samples <- max(hap.file.length, snp.file.length)
	
		if(snp.file.length < hap.file.length) {
	
			snp.repeats <-  rep(split.snp.file[[1]][length(split.snp.file[[1]])], (hap.file.length - snp.file.length))
			snp.repeats.append <- append(split.snp.file[[1]], snp.repeats)
			split.snp.file[[1]] <- snp.repeats.append
		}
	
	
		if(snp.file.length > hap.file.length) {
	
			hap.repeats <-  rep(split.hap.file[[1]][length(split.hap.file[[1]])], (snp.file.length - hap.file.length))
			hap.repeats.append <- append(split.hap.file[[1]], hap.repeats)
			split.hap.file[[1]] <- hap.repeats.append
		}
	
		gbl.var$hap.samples <- hap.file.length
	}
	
	split.color.list <- create.color.list(config.var)
	gbl.var$color.list <- split.color.list[[1]]
	
	split.symbol.fill.lists <- create.symbol.list(config.var, split.color.list, gbl.var)
	
	gbl.var$symbol.list <- split.symbol.fill.lists$split.symbol.list
	gbl.var$fill.list <- split.symbol.fill.lists$split.fill.list[[1]]
	

	#create sample labels
	#split SAMPLE.LABELS on comma, then take the top 5 or less entries and create a list
	if(!is.null(config.var$SAMPLE.LABELS)) {
		split.tmp.sample.list <- strsplit(config.var$SAMPLE.LABELS, ",")
	} else {
		split.tmp.sample.list <- rep("sample", samples)
	}
	
	if(length(split.tmp.sample.list[[1]]) < samples) {
		warning("SAMPLE.LABELS contains less labels than the number of samples. Labels must be separated by commas without spaces.\n")
	}
	
	split.tmp.sample.list.truncated <- substr(split.tmp.sample.list[[1]], 1, 12)
	
	#DEBUG STATEMENT
	#cat("samples ", samples, "\n")

	gbl.var$split.sample.labels <- list(head(split.tmp.sample.list.truncated, samples))
	
	grid.newpage()
	
	fix.var <- retrieve.data(split.snp.file, split.hap.file, config.var, gbl.var)
			
	config.var <- fix.var$config.var
	gbl.var <- fix.var$gbl.var

	if(config.var$DISP.LEGEND) {
		draw.legend(config.var, gbl.var)
	}

	#DEBUG STATEMENT 
	#cat("gbl.var ", is.null(gbl.var$geno.data), "\n")

	if(config.var$DISP.LDMAP) {
		draw.plot.ld.plot(config.var, gbl.var)
	}
	#END FUNCTION LIST

	#make sure the graphics device is closed and that the remaining device is the null device
	while(dev.cur()[[1]] != 1) {
		dev.off()
	}
		
	invisible(list(config.var, gbl.var))
}