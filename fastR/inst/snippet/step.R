###hop:3-9
step.model <- lm(HR - restHR ~ height * freq, step)
summary(step.model)
anova(step.model)
interaction.plot <- xyplot(HR - restHR ~ freq, data=step,
                           groups = height, type='a')
