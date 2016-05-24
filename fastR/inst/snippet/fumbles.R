m <- max(fumbles$week1)
table(factor(fumbles$week1,levels=0:m))
favstats(fumbles$week1)
