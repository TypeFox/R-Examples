binned.long2 <- cut(MammalLongevity$Longevity, breaks = c(0,10,20,30,40,50), right = FALSE)
tally( ~ binned.long2)      # no data frame given because it is not in a data frame

