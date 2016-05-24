library(dielectric)
library(ggplot2)

data(AgPalik)

AgPalik$set_span(200, 800)
raw <- dielectric2plot(AgPalik$raw())
silver <- AgPalik$predict(n=300, all.knots=TRUE, df=200)

d <- dielectric2plot(silver)

ggplot(d, aes(wavelength, value)) + geom_path() +
  facet_grid(variable~., scales="free") +
  geom_point(data=raw)



