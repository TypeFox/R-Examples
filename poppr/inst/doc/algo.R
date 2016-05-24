## ----echo=FALSE, message=FALSE-------------------------------------------
knitr::opts_knit$set(out.format = "latex")
thm <- knitr::knit_theme$get("acid")
knitr::knit_theme$set(thm)
knitr::opts_chunk$set(concordance=TRUE)
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
knitr::opts_chunk$set(out.width = '0.95\\linewidth', fig.align = "center", fig.show = 'asis')

## ------------------------------------------------------------------------
library("poppr")
dat.df <- data.frame(Genotype = c("1/1", "1/2", "2/3", "3/4", "4/4"))
dat <- as.genclone(df2genind(dat.df, sep = "/", ind.names = dat.df[[1]]))

## ------------------------------------------------------------------------
distances <- c("Nei", "Rogers", "Edwards", "Reynolds", "Prevosti")
dists <- lapply(distances, function(x){
  DISTFUN <- match.fun(paste(tolower(x), "dist", sep = "."))
  DISTFUN(dat)
})
names(dists) <- distances

# Adding Bruvo's distance at the end because we need to specify repeat length.
dists$Bruvo <- bruvo.dist(dat, replen = 1)
library("ape")
par(mfrow = c(2, 3))
x <- lapply(names(dists), function(x){ 
  plot(nj(dists[[x]]), main = x, type = "unrooted") 
  add.scale.bar(lcol = "red", length = 0.1)
})

