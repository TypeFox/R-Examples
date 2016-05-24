bwplot(Gender ~ TV, data = StudentSurvey)
dotPlot( ~ TV|Gender, layout = c(1, 2), width = 1, cex = 1, data = StudentSurvey)

