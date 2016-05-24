binned.long <- cut(MammalLongevity$Longevity, breaks = c(0,5,10,15,20,25,30,35,40))
tally( ~ binned.long)      # no data frame given because it is not in a data frame

