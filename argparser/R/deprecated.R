#' Create an argument parser.
#' 
#' This function is deprecated. Use \code{arg_parser} instead.
#'
#' @param description  description of the program
#' @param name         name of the program
#' @return a new \code{arg.parser} object
#' @export
arg.parser <- function(description, name=NULL) {
	warning("Use arg_parser() instead: arg.parser() is deprecated and will be removed in argparser-0.4.");
	arg_parser(description, name)
}

#' Add an argument to a parser.
#'
#' This function is deprecated. Use \code{add_argument} instead.
#' 
#' @param parser  an \code{arg.parser} object
#' @param arg     argument name (use no prefix for positional arguments,
#'                \code{--} or \code{-} prefix for optional arguments or flags)
#' @param help    help description for the argument
#' @param default default value for the argument [default: NA]
#' @param type    variable type of the argument (which can be inferred from 
#'                \code{default}), assumed to be \code{character} otherwise
#' @param flag    whether argument is a flag (and does not consume a value)
#'                [default: FALSE]
#' @param short   short-form for flags and positional arguments;
#'                short-forms can be assigned automatically based on the first
#'                character of the argument name, unless a conflict arises with
#'                an existing short-form; to avoid conflicts, add the argument 
#'                as early as possible
#' @return an \code{arg.parser} object with the argument added
#' @export
add.argument <- function(
	parser,
	arg, help,
	default=NULL, type=NULL, flag=NULL, short=NULL
) {
	warning("Use add_argument() instead: add.argument() is deprecated and will be removed in argparser-0.4.");
	add_argument(parser, arg, help, default, type, flag, short)
}

#' Parse arguments with a parser.
#'
#' This function is deprecated. Use \code{parse_args} instead.
#'
#' @param parser  an \code{arg.parser} object
#' @param argv    a character vector to parse (arguments and values should 
#'                already be split by whitespace)
#' @return a list with argument values
#' @export
parse.args <- function(parser, argv=commandArgs(trailingOnly=TRUE)) {
	warning("Use parse_args() instead: parse.args() is deprecated and will be removed in argparser-0.4.");
	parse_args(parser, argv)
}
