pulsePlot1 <- histogram(~pulse, data=littleSurvey)
pulsePlot2 <- histogram(~pulse, data=littleSurvey, subset=pulse>30) 
pulseSubset <- littleSurvey$pulse[ littleSurvey$pulse > 30 ]
mean(pulseSubset)
median(pulseSubset)
