#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2013 Robert Baron, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
#  
#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================

# CREATE BED ANNOTATION TRACK AND BEDGRAPH ANNOTATION TRACK
# READ IN bedplot.in OUTPUT FILE (bed.data.#) or ANY FILE CONTAINING MARKER, PHYSICAL POSITION AND SCORES

bedplot <- function(bed.data)
{

# Read in bed data file

if (file.exists(bed.data) == TRUE) { dat_bg <- read.table(bed.data, header=T) } else {

	warning(paste("File", bed.data, "does not", "exist!", sep=" "))
	return(FALSE)

}

dat_b <- dat_bg[dat_bg[,1]!="-", ]	

# Define the title of bedgraph table and bed table

t <- strsplit(bed.data, "\\.")[[1]][length(strsplit(bed.data, "\\.")[[1]])]

# Define the chromosome

if(strsplit(t, "0")[[1]][1] != "") chr <- t
if(strsplit(t, "0")[[1]][1] == "") chr <- strsplit(t, "0")[[1]][2]


n <- length(dat_bg) - 2

if (n == 1) {

	# Create bedgraph table
	
	w_bg <- dat_bg[!is.na(dat_bg[,3]),]

	# Define the length of bedgraph table

	n_bg <- length(w_bg[,1])

	if (n_bg != 0) {

		chr_bg <- rep(paste("chr", chr, sep=""), n_bg)

		end_bg <- start_bg <- w_bg[,2] - 1

		bg_line1 <- paste("browser position chr", chr, ":", min(start_bg), "-", max(start_bg), sep="") 
		bg_line2 <- paste("track type=bedGraph name=\"chr", chr, " results\" description=\"chr", chr, " results\" visibility=full color=200,100,0 altColor=0,100,200 priority=20", sep="")

		bg_t <- paste("BedGraph", t, sep=".")
		bg_tab <- cbind(chr_bg, start_bg, end_bg, w_bg[,3])

		write(bg_line1, file=bg_t, append=F)
		write(bg_line2, file=bg_t, append=T)
		write.table(bg_tab, file=bg_t, append=T, quote=F, row.names=F, col.names=F, sep="\t")

	}

	# Create bed table

	w_b <- dat_b[!is.na(dat_b[,3]),]

	# Define the length of bed table

	n_b <- length(w_b[,1])

	if (n_b != 0) {

		chr_b <- rep(paste("chr", chr, sep=""), n_b)
		start_b <- w_b[,2] - 1
		end_b <- start_b + 1
		strand <- rep("+", n_b)
		block <- rep("255,0,0", n_b)
	
		b_t <- paste("BED", t, sep=".")
		b_tab <- cbind(chr_b, start_b, end_b, as.character(w_b[,1]), w_b[,3]*100, strand, start_b, start_b, block)

		b_line1 <- paste("track name=\"chr", chr, " Stats\" description=\"Chr", chr, " stats\" visibility=3 useScore=1 itemRgb=\"On\" url=\"http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=$$\"", sep="")
		b_line2 <- paste("browser position chr", chr, ":", min(start_b), "-", max(start_b), sep="") 

		write(b_line1, file=b_t, append=F)
		write(b_line2, file=b_t, append=T)
		write.table(b_tab, file=b_t, append=T, quote=F, row.names=F, col.names=F, sep="\t")

	}

} else {

	for (i in 1:n) {

		# Create bedgraph table
	
		w_bg <- dat_bg[!is.na(dat_bg[,i+2]),]

		# Define the length of bedgraph table

		n_bg <- length(w_bg[,1])

		if (n_bg != 0) {

		chr_bg <- rep(paste("chr", chr, sep=""), n_bg)

		end_bg <- start_bg <- w_bg[,2] - 1

		bg_t <- paste("BedGraph", names(w_bg)[i+2], t, sep=".")
		bg_tab <- cbind(chr_bg, start_bg, end_bg, w_bg[,i+2])

		bg_line1 <- paste("browser position chr", chr, ":", min(start_bg), "-", max(start_bg), sep="") 
		bg_line2 <- paste("track type=bedGraph name=\"chr", chr, " results\" description=\"chr", chr, " results\" visibility=full color=200,100,0 altColor=0,100,200 priority=20", sep="")

		write(bg_line1, file=bg_t, append=F)
		write(bg_line2, file=bg_t, append=T)
		write.table(bg_tab, file=bg_t, append=T, quote=F, row.names=F, col.names=F, sep="\t")

		}

		# Create bed table

		w_b <- dat_b[!is.na(dat_b[,i+2]),]

		# Define the length of bed table

		n_b <- length(w_b[,1])

		if (n_b != 0) {

		chr_b <- rep(paste("chr", chr, sep=""), n_b)
		start_b <- w_b[,2] - 1
		end_b <- start_b + 1
		strand <- rep("+", n_b)
		block <- rep("255,0,0", n_b)
	
		b_t <- paste("BED", names(w_b)[i+2], t, sep=".")
		b_tab <- cbind(chr_b, start_b, end_b, as.character(w_b[,1]), w_b[,i+2]*100, strand, start_b, start_b, block)

		b_line1 <- paste("track name=\"chr", chr, " Stats\" description=\"Chr", chr, " stats\" visibility=3 useScore=1 itemRgb=\"On\" url=\"http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=$$\"", sep="")
		b_line2 <- paste("browser position chr", chr, ":", min(start_b), "-", max(start_b), sep="") 

		write(b_line1, file=b_t, append=F)
		write(b_line2, file=b_t, append=T)
		write.table(b_tab, file=b_t, append=T, quote=F, row.names=F, col.names=F, sep="\t")

		}		

	}

}

return(TRUE)

}
