# $Id: oddeven.R 1228 2007-11-30 18:05:43Z warnes $

# detect odd/even integers
odd <- function(x) x %% 2 == 1
even <- function(x) x %% 2 == 0 
