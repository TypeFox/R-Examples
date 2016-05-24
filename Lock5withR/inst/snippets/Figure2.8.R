histogram( ~ Pulse, type = "count", width = 5, data = StudentSurvey)
histogram( ~ Exercise, type = "count", width = 2, center = 2, 
          right = FALSE, ylim = c(0, 70), data = StudentSurvey) 
histogram( ~ Piercings, width = 1, data = StudentSurvey)

