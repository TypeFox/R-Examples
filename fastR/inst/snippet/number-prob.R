t <- table(littleSurvey$number); t
plot <- histogram(~number, littleSurvey, breaks=0.5 + 0:30)
max(t)
which(t == max(t))
which(t == min(t))
table(littleSurvey$number %% 2 == 0)
