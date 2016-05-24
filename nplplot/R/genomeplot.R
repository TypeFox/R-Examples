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

# CREATE TWO FORMAT GENOMEGRAPHS FILES CONTAINING EITHER CHROMOSOME BASE OR MARKER ID


genomeplot <- function(gg.data)
{

# Read in genome data file

if (file.exists(gg.data) == TRUE) { dat_chr <- read.table(gg.data, header=T) } else {

	warning(paste("File", gg.data, "does not", "exist!", sep=" "))
	return(FALSE)

}

dat_id <- dat_chr[dat_chr[,2]!="-", ]

n <- length(dat_chr) - 3

if (n == 1) {

	# Create chromosome table
	
	w_chr <- dat_chr[!is.na(dat_chr[,4]),]

	# Define the length of chromosome base table

	chr <- paste("chr", w_chr[,1], sep="")

	chr_tab <- data.frame(chr, w_chr[,3], w_chr[,4])
	names(chr_tab) <- c("Chromosome", "Position", names(w_chr)[4])

	write.table(chr_tab, file="GG.positions.all", append=F, quote=F, row.names=F, col.names=T, sep="\t")

	
	# Create id table

	w_id <- dat_id[!is.na(dat_id[,4]),]

	id_tab <- data.frame(as.character(w_id[,2]), w_id[,4])
	names(id_tab) <- c("Marker", names(w_id)[4])

	write.table(id_tab, file="GG.markers.all", append=F, quote=F, row.names=F, col.names=T, sep="\t")

} else {

	for (i in 1:n) {

		# Create chromosome table
	
		w_chr <- dat_chr[!is.na(dat_chr[,i+3]),]

		if (length(w_chr[,1]) != 0) {

			# Define the length of chromosome base table

			chr <- paste("chr", w_chr[,1], sep="")

			chr_t <- paste("GG.positions", names(w_chr)[i+3], "all", sep=".")
			chr_tab <- data.frame(chr, w_chr[,3], w_chr[,i+3])
			names(chr_tab) <- c("Chromosome", "Position", names(w_chr)[i+3])

			write.table(chr_tab, file=chr_t, append=F, quote=F, row.names=F, col.names=T, sep="\t")
		
		}
	
		# Create id table

		w_id <- dat_id[!is.na(dat_id[,i+3]),]

		if (length(w_id[,1]) != 0) {

			id_t <- paste("GG.markers", names(w_chr)[i+3], "all", sep=".")
			id_tab <- data.frame(as.character(w_id[,2]), w_id[,i+3])
			names(id_tab) <- c("Marker", names(w_id)[i+3])

			write.table(id_tab, file=id_t, append=F, quote=F, row.names=F, col.names=T, sep="\t")

		}

	}

}

return(TRUE)

}
