radviz2d <-
function (dataset, name = "") 
{
#    require(lattice)
    isfact = FALSE
    n = dim(dataset)[1]
    p = dim(dataset)[2]
    maintitle = paste("2D-Radviz for ", name)
    classes = dataset[, p]
    if (class(classes) == "factor") {
        isfact = TRUE
        classnames = levels(classes)
    }
    classnumbers = 1:length(unique(classes))
    classes = as.numeric(classes, drop = TRUE)
    varnames = colnames(dataset)
    dataset = as.matrix(mmnorm(dataset))
    dataset = dataset[, -p]
    sumrows = rowSums(dataset)
    columns = seq(0, (p - 2))
    angles = (2 * pi * columns)
    angles = angles/(p - 1)
    cosines = cos(angles)
    sines = sin(angles)
    proj.x = (dataset %*% cosines)
    proj.x = proj.x/sumrows
    proj.y = (dataset %*% sines)
    proj.y = proj.y/sumrows
    circledraw()
    anchors = p - 1
    def.par = par(font.lab = 2, font.sub = 2, cex = 0.6)
    par(xpd = TRUE)
    for (i in 1:anchors) {
        points(cosines[i], sines[i], pch = 19)
        if (cosines[i] > 0) 
            place = 4
        else if (cosines[i] < 0) 
            place = 2
        text(cosines[i], sines[i], varnames[i], pos = place, 
            font = 2)
    }
    for (j in 1:length(classnumbers)) {
        members = which(classes == classnumbers[j])
        for (i in members) {
            points(proj.x[i], proj.y[i], pch = (19), col = j + 
                1)
        }
    }
    if (!(isfact)) 
        classnames = paste("Class ", classnumbers)
    legend(x = -1.25, y = -1.25, legend = classnames, fill = 2:(length(classnumbers) + 
        1), bty = "n", xjust = 0, yjust = 0, horiz = TRUE)
    title(main = maintitle)
    par(xpd = FALSE)
    identify(proj.x,proj.y,labels=1:n)
    par(def.par)
}
