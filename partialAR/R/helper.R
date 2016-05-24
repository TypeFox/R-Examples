# helper.R -- A few miscellaneous helper functions
# Copyright (C) 2015 Matthew Clegg

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

ctext <- function (str, n) {
    # Centers the string str in a field of n blanks
    if (nchar(str) == n) return(str)
    if (nchar(str) > n) return(str, 1, n)
    
    nleft <- floor((n - nchar(str))/2)
    nright <- n - (nchar(str) + nleft)
    left_pad <- paste(rep(" ", nleft), collapse="")
    right_pad <- paste(rep(" ", nright), collapse="")
    pstr <- paste(left_pad, str, right_pad, sep="")
    pstr
}

printf <- function (...) { cat(sprintf(...)); }
println <- function (...) { cat(sprintf(...)); cat("\n"); }

vprintf <- function(fmt, X) {
  # Prints each element of a vector according to a specified format descriptor.
  # E.g., vprintf("%4.1f", c(1,2,3)) prints " 1.0 2.0 3.0"
  sapply (1:length(X), function(i) { cat(sprintf(fmt, X[i]))});
  NULL;
}

vprintln <- function(fmt, X) {
  # Same as vprintf but appends a newline character after the last item in
  # the vector has been printed.
  sapply (1:length(X), function(i) { cat(sprintf(fmt, X[i]))});
  cat ("\n");
  NULL;
}

quantile.table.from.samples <- function (colname, samples,
    quantiles = seq(0.01, 0.99, by=0.01)) {
    # Builds a table of quantile values from the data.frame samples
    
    sample_sizes <- sort(unique(samples$n))

    nq <- length(quantiles)
    ns <- length(sample_sizes)
    
    qtab <- matrix(NA, ncol=ns+1, nrow=nq+1)
    qtab[2:(nq+1), 1] <- quantiles
    qtab[1, 2:(ns+1)] <- sample_sizes
    
    for (i in 2:(ns+1)) {
        n <- sample_sizes[i-1]
        qtab[2:(nq+1), i] <- quantile(samples[samples$n == n, colname], quantiles)
    }

    qtab
}

quantile_table_interpolate <- function (qtab, sample_size, stat, stop.on.na=FALSE) {
	# On input, qtab is a dataframe of quantiles.  Each column corresponds to
	# a sample size, and each row corresponds to a quantile value.  The sample
	# sizes are given in the first row, and the quantiles are given in the
	# first column.  
	n <- nrow(qtab)
	i <- findInterval(sample_size, qtab[1,2:ncol(qtab)])+1
	if (i == 1) {
		parent_name <- as.character(sys.call(-1)[[1]])
		if (stop.on.na) {
			stop (parent_name," requires a minimum of ", qtab[1,2], " observations.")
		} else {
			warning (parent_name, " requires a minimum of ", qtab[1,2], " observations.")
			return(NA)
		}
	}
	y1 <- approx(qtab[2:n, i], qtab[2:n, 1], stat, rule=2)$y
	if (i < ncol(qtab)) {
		y2 <- approx(qtab[2:n, i+1], qtab[2:n, 1], stat, rule=2)$y
		n1 <- qtab[1,i]
		n2 <- qtab[1,i+1]
		y <- y1 * (n2 - sample_size) / (n2 - n1) + y2 * (sample_size - n1)/(n2 - n1)
	} else {
		y <- y1
	}
	y
}

