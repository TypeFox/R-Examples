library(dielectric)
library(ggplot2)

data(AuJC)


AuJC$set_span(300, 800)
raw <- dielectric2plot(AuJC$raw())
gold <- AuJC$predict(n=300, all.knots=TRUE)

d <- dielectric2plot(gold)

ggplot(d, aes(wavelength, value)) + geom_path() +
  facet_grid(variable~., scales="free") +
  geom_point(data=raw)



