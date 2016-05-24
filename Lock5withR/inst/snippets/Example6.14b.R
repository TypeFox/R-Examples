favstats(~Time, data = CommuteAtlanta)
confint(t.test(~Time, conf.level = 0.99, data = CommuteAtlanta))
confint(t.test(~Time, conf.level = 0.95, data = CommuteAtlanta))

