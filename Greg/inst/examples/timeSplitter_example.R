test_data <- data.frame(
  id = 1:4,
  time = c(4, 3.5, 1, 5),
  event = c("alive", "censored", "dead", "dead"),
  age = c(62.2, 55.3, 73.7, 46.3),
  date = as.Date(
    c("2003-01-01", 
      "2010-04-01", 
      "2013-09-20",
      "2002-02-23"))
)
timeSplitter(test_data, .5, 
             time_var = "time",
             time_related_vars = c("age", "date"),
             event_var = "event")
