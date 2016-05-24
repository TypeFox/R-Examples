#!/usr/bin/env Rscript

library(argparser, quietly=TRUE)
    
# Create a parser
p <- arg_parser("Calculate the sum of a set of numbers")

# Add command line arguments
p <- add_argument(p, "--numbers", help="numbers to add",
	nargs=Inf, type="numeric")
p <- add_argument(p, "--subset",
	help="start and end (1-based) indices of subset of numbers to add",
	default=c(1, Inf));

# Parse the command line arguments
argv <- parse_args(p)

# Do work based on the passed arguments
idx <- argv$subset;
if (is.infinite(idx[2])) {
	idx[2] <- length(argv$numbers);
}
cat( sum(argv$numbers[idx[1]:idx[2]]), "\n")

