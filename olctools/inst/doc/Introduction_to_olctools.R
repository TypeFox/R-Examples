## ---- eval=FALSE---------------------------------------------------------
#  library(olctools)
#  # Encode the location of Mbgathi Road, in Nairobi, Kenya
#  encode_olc(-1.314063, 36.79881, 10)
#  
#  # [1] "6GCRMQPX+9G"

## ---- eval=FALSE---------------------------------------------------------
#  shorten_olc("6GCRMQPX+9G", -1.314063, 36.79881)
#  
#  # [1] "+9G"

## ---- eval=FALSE---------------------------------------------------------
#  str(decode_olc(c("6GCRMQPX+9G", "7FG49QCJ+2VX")))
#  
#  # 'data.frame':    2 obs. of  7 variables:
#  #  $ latitude_low    : num  -1.31 20.37
#  #  $ longitude_low   : num  36.8 2.78
#  #  $ latitude_center : num  -1.31 20.37
#  #  $ longitude_center: num  36.8 2.78
#  #  $ latitude_high   : num  -1.31 20.37
#  #  $ longitude_high  : num  36.8 2.78
#  #  $ code_lengths    : int  10 11

## ---- eval=FALSE---------------------------------------------------------
#  recover_olc("9G8F+6X", 47.4, 8.6)
#  
#  # [1] "8FVC9G8F+6X"

