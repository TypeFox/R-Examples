library(adegraphics)
pdf("slogo.pdf")

## ex1
data(ggtortoises, package = "ade4")
ico <- ggtortoises$ico[as.character(ggtortoises$pop$carap)]
g1 <- s.logo(ggtortoises$pop, ico, pori.incl = FALSE)
g2 <- s.label(ggtortoises$pop, add = TRUE, plabels = list(boxes = list(alpha = 0.4, border = "transparent")))

## ex2
data(capitales, package = "ade4")
index <- unlist(lapply(1:15, function(i) which(names(capitales$logo) == tolower(rownames(capitales$xy)[i]))))
g3 <- s.logo(capitales$xy, capitales$logo[index])

x <- c(0, max(capitales$area$x))
y <- c(0, max(capitales$area$y))
#g4 <- s.image(cbind(x, y), z = c(1, 2), outsideLimits = capitales$area, grid = 500, regions = list(col = "yellow", alpha = 0.9))
#s.logo(capitales$xy, capitales$logo[index], add = TRUE)
