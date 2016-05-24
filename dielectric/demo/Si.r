library(dielectric)
library(ggplot2)

data(aSi)

aSi$set_span(300, 800)
raw <- dielectric2plot(aSi$raw())
silicon <- aSi$predict(n=300, all.knots=TRUE)

d <- dielectric2plot(silicon)

ggplot(d, aes(wavelength, value)) + geom_path() +
  facet_grid(variable~., scales="free") +
  geom_point(data=raw)



