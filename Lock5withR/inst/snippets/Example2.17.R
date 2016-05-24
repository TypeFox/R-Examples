histogram( ~ Pulse, fit = "normal", data = StudentSurvey)
mean <- mean( ~ Pulse, data = StudentSurvey); mean
sd <- sd( ~ Pulse, data = StudentSurvey); sd
mean - 2*sd
mean + 2*sd

