`getHH` <-
function (YYYYMMDDHH)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
# require("chron")

 YYYYMMDDHH <- sapply(YYYYMMDDHH, as.character)
 l <- sapply(YYYYMMDDHH, nchar)
 if (any(I <- (l > 10 | l == 9 | l < 8))) stop("invalid date string")
 HH <- sapply(YYYYMMDDHH, 
   function(s) if (nchar(s) == 10) substring(s,9,10) else"00")
 unique(HH)
}

