################################################################################
# checks.R
# (c)2012 Andrew Redd 
# This is file part of the lint R package, a code style check package for R.
# 
# This file contains the examples of bad code for test validation.
# Checks validate by line number starting at 11. do NOT alter the line order!
################################################################################

very.long.line.80 <- function(){print("hello world for a very long line, that keeps going.")}
very.long.line.100 <- function(){print("hello world for a very long line, that keeps going past 100.")}
	# literal tab
   # odd indent
do.call(
   with.odd.indent)
c(1,2)  # no space after comma.
# ignore no space after comma,in comment.
c(1,
  2)  # end of lines do not trigger check
  
{  # spacing.twobeforecomments
lm# needs two spaces
lm # needs two spaces
lm  # This one OK
}
