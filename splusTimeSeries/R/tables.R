axis.break.table <-
  list(
       BreakSpan = timeSpan(
         julian = c(18262, 2191, 1096, 365, 60, 7, 1, 0, 0, 0),
         ms = 60*60*1000*c(0, 0, 0, 0, 0, 0, 0, 12, 6, 6)),
       SampleTime = timeSpan(
         julian = c(365, 90, 30, 7, 1, 0, 0, 0, 0, 0),
         ms = 60*1000*c(0, 0, 0, 0, 0, 3*60, 60, 30, 15, 5)),
       Align.by = c("years", "quarters", "months", "weeks", "days",
         "hours", "hours", "minutes", "minutes", "minutes"),
       Align.k.by = c(1, 1, 1, 1, 1, 3, 1, 30, 15, 15)
       )

axis.label.table <-
  list(
       Units = c("time.of.day", "days", "days", "days", "days", "weeks",
         "weeks", "weeks", "months", "months", "quarters", "quarters", "years"),
       Number = c(1, 1, 2, 3, 4, 3, 2, 1, 1, 2, 1, 2, 1),
       Label = c("time.of.day", "days", "days", "month.day", "date", "date",
         "days", "month.day", "months", "month.year", "quarters",
         "quarter.year", "years"),
       Outer.Label = c("date", "date", "month.year", "years", "", "",
         "month.year", "years", "years", "", "years", "", ""),
       Outer.By = c("days", "weeks", "months", "years", "", "", "months",
         "years", "years", "", "years", "", ""))

axis.tick.table <-
  list(
       table = c(rep(1, 12), rep(2, 5), rep(3, 4), rep(4, 3), rep(5, 4),
         rep(6, 6)),
       small.k.by = c(1, 1, 2, 5, 10, 25, 50, 100, 250, 500,
         1, 5, 15, 1, 1, 5, 5, 15, 15, 1, 3, 3, 6, rep(1, 8), 5, 25, 100),
       small.by = c(rep("milliseconds", 10), rep("seconds", 3),
         rep("minutes", 6), rep("hours", 4), rep("days", 1),
         rep("weeks", 2), rep("months", 2), rep("quarters", 2),
         rep("years", 4)),
       medium.k.by = c(5, 5, 10, 25, 50, 100, 250, 500, 1, 2, 5, 15,
         1, 5, 5, 15, 15, 1, 1, 6, 12, 12, rep(1, 8), 5, 10, 500, 500),
       medium.by = c(rep("milliseconds", 8), rep("seconds", 4),
         rep("minutes", 5), rep("hours", 5), rep("days", 1),
         rep("weeks", 1), rep("months", 2), rep("quarters", 2),
         rep("years", 6)),
       big.k.by = c(10, 25, 100, 250, 500, 1, 5, 5, 10, 10, 15, 1, 5, 15, 1,
         1, 3, 3, 6, rep(1, 8), 5, 5, 5, 10, 25, 1000, 5000),
       big.by = c(rep("milliseconds", 5), rep("seconds", 6),
         rep("minutes", 3), rep("hours", 5), rep("days", 2),
         rep("weeks", 2), "months", "quarters", rep("years", 9)),
       span = timeSpan(
         julian = c(rep(0,23), 1, 7, 7, 30, 30, 91, 91, 365, 1826, 9131, 36525),
         ms = c(1, 1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 5*1000, 15*1000,
           60*1000, 60*1000, 5*60*1000, 5*60*1000, 15*60*1000, 15*60*1000,
           60*60*1000, 3*60*60*1000, 3*60*60*1000, 6*60*60*1000, 0, 0, 0,
           (10*60 + 30)*60*1000, (10*60 + 30)*60*1000, (7*60 + 30)*60*1000, 
           (7*60 + 30)*60*1000, 6*60*60*1000, 6*60*60*1000, 6*60*60*1000, 0)
         ),
       span.frac = c(0.2, rep(1, 12), 0.5, 1, 0.5, 1, 0.5, rep(1, 6),
         0.5, 1, 0.3, 1, 0.5, rep(1, 5)))
