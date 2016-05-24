library(dyn)
data(airquality)
start.date <- as.Date(paste(c(1973, airquality[1,5:6]), collapse = "-")) 
airquality.zoo <- zooreg(data.matrix(airquality), start = start.date)

ozone.lm <- dyn$lm(Ozone ~ ., airquality.zoo, na.action = na.omit)
ozone.lm
head(fitted(ozone.lm))
head(resid(ozone.lm))

ozone2.lm <- dyn$lm(Ozone ~ ., airquality.zoo, na.action = na.exclude)
ozone2.lm
head(fitted(ozone.lm))
head(resid(ozone.lm))

end. <- time(airquality.zoo)[140]

# subset
ozone3.lm <- dyn$lm(Ozone[time(Ozone) < end.] ~ lag(Ozone, -1), airquality.zoo)
ozone3.lm

# alternative subset
ozone4.lm <- dyn$lm(window(Ozone, end = end.) ~ lag(Ozone, -1), airquality.zoo)
ozone4.lm

