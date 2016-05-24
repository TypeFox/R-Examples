traffic

reshape(traffic[,-2], idvar="year",ids=row.names(traffic),
        times=names(traffic)[3:6],timevar="state",
        varying=list(names(traffic)[3:6]),
        v.names="deathRate",
        direction="long") -> longTraffic
head(longTraffic)
