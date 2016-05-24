head(CommuteStLouis)
favstats(~Time, data = CommuteStLouis)
favstats(~Time, data = CommuteAtlanta)
bwplot(~Time, xlim = c(0, 200), data = CommuteAtlanta) # to check for normality
bwplot(~Time, xlim = c(0, 200), data = CommuteStLouis) # to check for normality

