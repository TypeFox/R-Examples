airp.aov <- aov(pollution~location,airpollution)
TukeyHSD(airp.aov)
