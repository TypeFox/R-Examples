head(StudentSurvey, 3)
favstats( ~ Pulse, data = StudentSurvey)
favstats(Pulse ~ Award, data = StudentSurvey)
anova(lm(Pulse ~ Award, StudentSurvey))

