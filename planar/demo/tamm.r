# ---- setup -----
library(planar)
opts_chunk$set(fig.path="tamm/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))

# prettier palette
palette(palette_tamm)

scale_colour_discrete <- function(...) 
  scale_colour_brewer(..., palette="Set1")

scale_fill_discrete <- function(...) 
  scale_fill_manual(..., values=palette())

# ---- structure -----

tamm <- tamm_stack_ir(pairs=13, dm=100)
(p <- autoplot(tamm))

# ---- ff -----
ff <- simulate_ff(s=tamm, wavelength=seq(600, 1200))
head(ff)
mff <- melt(ff, meas=c("R","T","A"))
ggplot(mff, aes(wavelength, value, colour=variable))+
  geom_line() +
  scale_y_continuous(lim=c(0,1), expand=c(0,0)) +
  labs(x = "wavelength /nm", y="")

optimum <- subset(ff, A == max(A))
optimum

# ---- nf -----
nf <- simulate_nf(s=tamm, wavelength=optimum$wavelength)
head(nf)

p + geom_line(aes(x, I), data=nf) +
  scale_y_continuous(expression("|E|"^2), expand=c(0,0))

