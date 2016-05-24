confint(t.test(Wetsuits$Wetsuit, Wetsuits$NoWetsuit, paired = TRUE))
confint(t.test( ~ (Wetsuit - NoWetsuit), data = Wetsuits))

