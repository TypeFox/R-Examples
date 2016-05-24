# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzelez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


bin.color <-
function(mat, threshold, breaks, col, symkey) 
{
    if (isTRUE(symkey)) {
        max.mat = max(abs(mat))
        min.mat = -max.mat
    }
    else {
        max.mat = max(mat)
        min.mat = min(mat)
    }
	
	    if (missing(breaks) || is.null(breaks)) {
        if (class(col) == "function") breaks = 32
            else breaks = length(col) 
    }

    if (length(breaks) == 1) {
		if (isTRUE(symkey)) {
				if ((breaks/2) - trunc(breaks/2) != 0) 
					stop("'breaks' must be a even number if 'symkey = TRUE'", call. = FALSE)
			
			if (threshold == 0) breaks = c(seq(min.mat, max.mat, length = breaks + 1))
			else {			
				nb = breaks/2
				breaks = c(seq(min.mat, -threshold, length = nb + 1), 0, 
						   seq(threshold, max.mat, length = nb + 1))					   
				id = which(breaks == 0)
				breaks = breaks[-c(id - 1, id + 1)]			
			}
		}
		else { 
			breaks = breaks + 1
			
			if ((min.mat < -threshold) & (max.mat < threshold))
				breaks = seq(min.mat, -threshold, length = breaks)
				
			if ((min.mat > -threshold) & (max.mat > threshold))
				breaks = seq(threshold, max.mat, length = breaks)
				
			if ((min.mat < -threshold) & (max.mat > threshold)) {
				if (threshold == 0) breaks = c(seq(min.mat, max.mat, length = breaks))
				else {
					long = max.mat - min.mat - 2*threshold
					bin = long/breaks
					breaks = seq(threshold, -min.mat, by = bin)
					o = order(breaks, decreasing = TRUE)
					breaks = c(-breaks[o], 0, seq(threshold, max.mat, by = bin))
					id = which(breaks == 0)				
					breaks = breaks[-c(id - 1, id + 1)]
				}
			}
		}
    }

    ncol = length(breaks) - 1

    if (class(col) == "function") 
        col = col(ncol)
		
	if (length(breaks) != length(col) + 1) 
        stop("must have one more break than colour", call. = FALSE)	

    min.breaks = min(breaks)
    max.breaks = max(breaks)

    mat[mat < min.breaks] = min.breaks
    mat[mat > max.breaks] = max.breaks
    
    bin <- .bincode(as.double(mat), as.double(breaks), TRUE, TRUE)

    return(invisible(list(bin = bin, col = col, breaks = breaks, lim = c(min.mat, max.mat))))		
}
