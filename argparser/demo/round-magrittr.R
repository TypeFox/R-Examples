#!/usr/bin/env Rscript

library(magrittr, quietly=TRUE)
library(argparser, quietly=TRUE)
    
# Create a parser and add arguments, and parse the command line arguments
# All steps are chained together, courtesy of magrittr
p <- arg_parser("Round a floating point number") %>%
	add_argument("number", help="number to round", type="numeric") %>%
	add_argument("--digits", help="number of decimal places", default=0)

# Parse the command line arguments
# This step is kept separate to simply error message
argv <- parse_args(p)

# Do work based on the passed arguments
cat( round(argv$number, argv$digits), "\n")

