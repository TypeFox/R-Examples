x.bar <- mean(~Time, data = CommuteAtlanta); x.bar
SE <- sd( ~ mean, data = Bootstrap ); SE        # standard error
MoE <- 2 * SE; MoE                              # margin of error for 95% CI
x.bar - MoE                                     # lower limit of 95% CI
x.bar + MoE                                     # upper limit of 95% CI

