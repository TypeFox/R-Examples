library(parfm)
data(culling)
head(culling)
culling <- culling[culling$Time > 0,]
culling$TimeMonths <- culling$Time * 12 / 365.25

set.seed(1)
culling <- culling[sample(1:nrow(culling), 150),]
parfm(Surv(TimeMonths, Status)~LogSCC, cluster="Herd", data=culling,
               dist="exponential", frailty="gamma")
