# Verified 1.3.18
# Version 4.0
panomapa <- function(collection, main, axis = TRUE,
                   xlab = "Long",
                   ylab = "Lat",
                   lab.col = 'black',
                   bg = NA,
                   map.bg = NA,
                   map.col = 'black',
                   col.ramp = c("Green3", "darkorange1","red"),
                   arrow.cex = 4.5,
                   arrow.plot = TRUE,
                   pt.col = rgb(0, 0, 0, 0.75),
                   pt.cex = 4.5,
                   pt.pch = 21,
                   leg.pt.bg = pt.bg,
                   leg.bg = NA,
                   leg.title = "Lengevity\n(years)",# "Longevidad\n(años)",
                   leg.offset = c(0,0),
                   leg.y.intersp = 1.75) {


colf = function(x) {
  colorRamp(col.ramp)(x)
}

print_ramp <- function (ColorLevels, pal, mar = c(1,2.5,2.5,2), xlab="", ylab="", ...) {
  par(mar=mar)
  
  image(1, ColorLevels,
        matrix(data = ColorLevels, ncol = length(ColorLevels), nrow = 1),
        col = pal,
        xaxt = "n", ylab = ylab, xlab = xlab, ...)
  title(...)
}

# Function start ####

if ( length(names(collection$catalog)) != 0 ) {
        cat = collection$catalog; catalogo = list(); catalogo[[1]] = cat; rm("cat")
        dat = collection$data; datos = list(); datos[[1]] = dat; rm("dat")
} else {
        catalogo = collection$catalog
        datos = collection$data
}
map.abb = unique(unlist(lapply( catalogo, function(x){x$State} )))

if ( length(catalogo) != length(datos) ) {
  stop("In collection: catalogo and datos lengths differ.")
}
n = length(catalogo)
pos = plyr::ldply(catalogo, function(x) {
  c(x$Longitude, x$Latitude)
})
D = matrix(c(range(pos[,1]), range(pos[,2])), ncol=2)
dpa = list()
for (k in 1:n) {
  dpa[[k]] = list(qty.disp = length(catalogo[[k]]$Avble.yrs), 
                  m = length(datos[[k]]), frac.na = sum(is.na(datos[[k]]))/length(datos[[k]]), 
                  frac.ag = sum(datos[[k]] < 0, na.rm = T)/length(datos[[k]]))
}
pta = plyr::ldply(dpa, function(x) {
  c(x$qty.disp, x$m, x$frac.na, x$frac.ag)
})
m.disp = max(pta[, 1])
pta[, 1] = pta[, 1]/m.disp
pta[pta[, 1] <= 0.1, 1] = 0.1

pt.bg = rgb(colf(pta[, 3])/255, alpha = 0.75)
leg.pt.bg = rgb(colf(rev(c(0.1, 0.25, 0.5, 0.75, 1)))/255, 
                alpha = 0.75)
par.save <- par(no.readonly = TRUE)
layout(mat = matrix(c(1,1,1,2,3,4), ncol=2), widths = c(4,1), heights = c(2,2,1))

ptatruescale = par()$fin[1] * 0.66

if (!is.na(map.abb)) {
        ESS <- get.shape.state(map.abb)
        SHP.range = matrix(ncol = 4, nrow = length(ESS))
        for (i in 1:length(ESS)) {
                d = slot(ESS, "polygons")[[i]]
                SHP.sub = matrix(ncol = 4, nrow = length(slot(d,"Polygons")))
                for (j in 1:length(slot(d, "Polygons"))) {
                        d.sub = slot(d, "Polygons")[[j]]
                        d.sub = slot(d.sub, "coords")
                        SHP.sub[j, 1:2] = range(d.sub[, 1])
                        SHP.sub[j, 3:4] = range(d.sub[, 2])
                }
                d = matrix(apply(SHP.sub, 2, range), ncol = 4)
                SHP.range[i, 1:2] = diag(d[1:2, 1:2])
                SHP.range[i, 3:4] = diag(d[1:2, 3:4])
        }
        d = matrix(apply(SHP.range, 2, range), ncol = 2)
        D = rbind(d, D) # max betwen points and shape border
}

plot(axes = F, asp = 1, bty = "n", type = "n", range(D[,1]), range(D[, 2]), ylab = ylab, xlab = xlab)
if (!is.na(map.abb)) {
        plot(add = T, axes = F, ESS, bg = map.bg, border = map.col, asp = 1)
}

points(pos, cex = ptatruescale * pt.cex * pta[, 1], bg = pt.bg, pch = pt.pch, 
       col = pt.col) 

if (axis == T) {
        axis(1, col = lab.col, col.axis = lab.col)
        axis(2, col = lab.col, col.axis = lab.col)
}

if (missing(main)) {
        if ( length(catalogo) == 1) {
                main = "Station longevity"
        } else {
                main = "Stations longevity"
        }
        if ( ! is.na(map.abb) ) {
                estados.venezuela <- get.shape.state()
                main = paste(main, "for", paste(estados.venezuela[map.abb, "shape.name"], collapse = ", "))
        }
}
title(main = main, col.main = lab.col,cex.main=2.5)
long = round(c(0.1, 0.25, 0.5, 0.75, 1) * m.disp, 0)
long = apply(cbind(c("<", "<", "<", "<", "<"), long), 1, paste0, collapse = "")
par(mar = c(0.5,0.5,0,0.5) + 0.1, mai=c(0,0,1,0))
plot(c(-1,1), c(-1,6), typ='n', asp=1, axes=F, xlab=NA, ylab=NA)
legend(x = -1, y = 5.9,  
       legend = long, 
       pt.cex = pt.cex * ptatruescale * c(0.2,0.25, 0.5, 0.75, 1), 
       pch = 21, bg = leg.bg, pt.bg = NA,
       cex = 1.25, bty = "n", text.col = lab.col, 
       y.intersp = leg.y.intersp, )
title(main = leg.title, cex.main = 1.45, font.main = 2)
leg.lvl = seq(0, 100, by=5)
leg.col = rgb(colf(rev(leg.lvl/100))/255, alpha = 0.75)
print_ramp(leg.lvl, leg.col, main="Data %",mar = c(1,5,7.5,3.5))

par(mar = rep(0.5,4) + 0.1)
plot(c(-1,1), c(-1,5), typ = 'n', asp = 1, axes = F, xlab = NA, ylab = NA)
if (arrow.plot){ plotArrow(cex = arrow.cex) }
par(par.save)
}