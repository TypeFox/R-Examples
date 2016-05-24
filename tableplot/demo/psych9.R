# psych9.R: Factor analysis of 9 psychological abilities from datasets::Harman74.cor
# 
# Author: Michael Friendly
###############################################################################

v9 <-c(1,2,4,6,7,9,10,12,13)
psych9<-Harman74.cor

# factor analysis, with varimax rotation
fa3 <- factanal(covmat=psych9,factors=3)

unique <- fa3$uniquenesses
communalities <- 1- fa3$uniquenesses

vnames <- rownames(fa3$loadings)
snames <- abbreviate(vnames,8)

# assign factor names according to interpretation
fnames <- c("Verbal", "Visual", "Speed")

loadings <- fa3$loadings

# override print method for loadings
class(loadings)<-"matrix"
colnames(loadings) <- fnames
rownames(loadings) <- snames

# display horizontally, rounded
values <- round(100 * t(loadings));

## choose colors and symbols to highlight the assumed interpretation

specs <- make.specs(
		shape=c(0, 0, 0, 2),           #  circles and squares=non-target loadings
		shape.col.neg="black",
		cell.fill=c("red","blue","green", "grey40"),
		back.fill="white",
		scale.max=100,
		label=1, label.size=1, 
)

types <- matrix(
        c(rep(c(4, 2, 4), each=3),
		  rep(c(1, 4, 4), each=3),
          rep(c(4, 4, 3), each=3)), 3, 9, byrow=TRUE)

tableplot(values, assign.sets=types,
		cell.specs=specs,
		title="Tableplot of varimax rotated loadings",
		side.rot=90,
		top.space=20, left.space=8
)

# add a row showing unique variances

cvalues <- rbind(values, round(100*unique))
rownames(cvalues) <- c(fnames,"Unique")
cspecs <- make.specs(
		shape=c(0, 0, 0, 2, 2),           #  circles and squares
		shape.col.neg="black",
		cell.fill=c("red","blue","green", "grey40", "white"),
		back.fill=c(rep("white",4), "grey"),
		scale.max=100,
		label=1, label.size=1, 
)
ctypes <- rbind(types, 5)

tableplot(cvalues, ctypes,
		title="Tableplot of varimax rotated loadings",
		cell.specs=cspecs,
		side.rot=90,
		h.parts=c(3,1), gap=4,
		top.space=20, left.space=8,
)

## try promax rotation
ploadings <- promax(loadings(fa3),m=3)$loadings
class(ploadings)<-"matrix"
colnames(ploadings) <- fnames
rownames(ploadings) <- snames
ploadings

values <- round(100 * t(ploadings));

tableplot(values, types,
		cell.specs=specs,
		title="Tableplot of promax rotated loadings",
		side.rot=90,
		top.space=20, left.space=8,
)


