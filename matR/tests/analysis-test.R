library(matR)
N <- 1:4
List <- mget (paste0 ("xx", N), inherits=TRUE)

#-----------------------------------------------------------------------------------------
#  OK FOR CRAN
#-----------------------------------------------------------------------------------------

for (xx in List) {
#-----------------------------------------------------------------------------------------
#  distx()		...biom method
#-----------------------------------------------------------------------------------------
	uu <- as.matrix (xx, TRUE) [,1]
	vv <- as.matrix (xx, TRUE) [1,]
	distx (xx)														# distance between columns
	distx (xx, bycol=FALSE)											# distance between rows
	distx (xx, method="bray-curtis")								# alt measure
	distx (xx, method="bray-curtis", bycol=FALSE)
	distx (xx, groups=1:ncol(xx) %% 4)								# mean pairwise distance between groups
	distx (xx, groups=1:nrow(xx) %% 4, bycol=FALSE)					# row groups
	distx (xx, p=uu)												# from each col to a given vector
	distx (xx, p=vv, bycol=FALSE)									# from each row
	distx (xx, p=uu, groups=1:ncol(xx) %% 4)						# from each group to given vector
	distx (xx, p=vv, groups=1:nrow(xx) %% 4, bycol=FALSE)			# row groups
	}

for (xx in List) {
#-----------------------------------------------------------------------------------------
#  rowstats()	...biom method
#-----------------------------------------------------------------------------------------
	str (rowstats (xx, groups=seq(along=colnames(xx)) %% 2, test="Kr"))
	str (rowstats (xx, groups=seq(along=colnames(xx)) %% 3, test="Kr"))
	str (rowstats (xx, groups=seq(along=colnames(xx)) %% 2, test="t-test-un"))
	str (rowstats (xx, groups=seq(along=colnames(xx)) %% 2, test="Mann"))		# gives warning re. ties
	str (rowstats (xx, groups=seq(along=colnames(xx)) %% 2, test="AN"))
	str (rowstats (xx, groups=seq(along=colnames(xx)) %% 3, test="AN"))
	if (ncol(xx) %% 2 != 1) {
		str (rowstats (xx, groups=seq(along=colnames(xx)) %% 2, test="t-test-p"))
		str (rowstats (xx, groups=seq(along=colnames(xx)) %% 2, test="Wilc"))
		}
	}

for (xx in List) {
#-----------------------------------------------------------------------------------------
#  transform()				
#-----------------------------------------------------------------------------------------
	transform (xx, t_NA2Zero)
	transform (xx, t_NA2Zero, t_Threshold)
	transform (xx, t_NA2Zero, t_Threshold = list(entry=3))
	transform (xx, t_NA2Zero, t_Threshold = list(row=6))
	transform (xx, t_NA2Zero, t_Threshold = list(row=6,col=9))
	transform (xx, t_NA2Zero, t_Threshold = list(entry=5))
	transform (xx, t_NA2Zero, t_Threshold = list(entry=5), t_Log)
	transform (xx, t_NA2Zero, t_Threshold = list(entry=5), t_Log, t_ColCenter)
	}

for (xx in List) {
#-----------------------------------------------------------------------------------------
#  boxplot()
#-----------------------------------------------------------------------------------------
	xx.normed <- transform (xx, t_Log)
	boxplot(xx)
	boxplot(
		xx, 
		xx,
		main="so good they named it twice")
	boxplot(
		xx, 
		xx.normed)
	boxplot( 
		xx, 
		xx.normed,
		columns=2:4)
	boxplot(
		xx, 
		xx.normed,
		x.main="raw",
		y.main="log")
	boxplot(
		xx,
		xx.normed,
		x.main="raw",
		y.main="log",
		cex.main=2,
		x.names="$$project.id",
		x.cex.axis=1.5,
		y.names="$$metagenome.id",
		y.cex.axis=0.75)
	}
xx <- xx1 ; xx.normed <- transform (xx, t_Log)
boxplot(
	xx,
	xx.normed,
	map=c(
		col="host_common_name"))
boxplot(
	xx, 
	xx.normed,
	x.main="raw",
	y.main="log",
	map=c(
		col="host_common_name"))
boxplot(
	xx.normed,
	xx.normed,
	x.main="log",
	y.main="log",
	map=c(
		x.col="host_common_name",
		y.col="samp_stor"),
	y.col=c(
		"-80"="salmon",
		"NA"="orange"))

#-----------------------------------------------------------------------------------------
#  princomp()				
#-----------------------------------------------------------------------------------------
princomp (xx1, method="euclidean")
princomp (xx1, method="bray-curtis")

princomp (xx1, dim=1)											# single PC
princomp (xx1, dim=2)
princomp (xx1, dim=c(1,2))										# two PCs
princomp (xx1, dim=c(2,3))
princomp (xx1, dim=c(1,2,3))									# same three PCs
princomp (xx1, dim=c(1,3,2))									# from different perspectives
princomp (xx1, dim=c(2,1,3))
princomp (xx1, dim=c(2,3,4))									# different three PCs

princomp (xx1, labels = "")										# labeling variations (color, size, metadata)
princomp (xx1, labels = LETTERS [1:7])
princomp (xx1, labels = LETTERS [1:7], label.col = "blue")
princomp (xx1, labels = "$$host_common_name")
princomp(
	xx1, 
	labels="$$pubmed_id",
	label.col="blue",
	label.cex=0.5)
# princomp(
# 	xx1,
# 	labels="$$pubmed_id", 
# 	map=c(
# 		label.col="host_common_name"))
# princomp(
# 	xx1,
# 	labels="$$pubmed_id", 
# 	map=c(
# 		label.col="host_common_name"),
# 	label.col=c(
# 		"cow"="blue",
# 		"striped bass"="brown",
# 		"Mouse"="brown"))
princomp (xx3, dim=3, labels="", map=c(col="biome"))

princomp(														# plotting character variations
	xx1,
	col="blue")
princomp(
	xx1,
	col=c("blue","blue","blue","red","red","red","red"))
princomp(
	xx1,
	pch=17)
princomp(
	xx1,
	pch=15:21)
princomp(
	xx1,
	cex=2)
princomp(
	xx1,
	cex=seq(1,2,len=7))

princomp(
	xx1, 
	map=c(												# automap one par variable to metadata
		pch="samp_store_temp"))
princomp(
	xx1, 
	map=c(												# automap two
		col="host_common_name",
		pch="samp_store_temp"))
princomp(
	xx1, 
	map=c(												# automap three
		col="host_common_name",
		pch="samp_store_temp",
		cex="material"))
princomp (xx1,
	map=c(												# explicitly map one (of two)
		col="host_common_name",
		pch="samp_store_temp"),
	col=c(
		Mouse="brown",
		cow="red",
		"striped bass"="blue"))
princomp (xx1,
	map=c(												# explicitly map both
		col="host_common_name",
		pch="samp_store_temp"),
	col=c(
		Mouse="brown",
		cow="red",
		"striped bass"="blue"),
	pch=c(
		"-80"="+",
		"NA"="x"),
	cex=2)
princomp(
	xx1,
	dim=1:2,
	map=c(												# give explicit but incomplete map
		cex="host_common_name"),
	cex=c(
		cow=2.5),
	labels="")

zz <- princomp (xx1)									# reuse computation
princomp(
	xx1, 
	main="title added with\nno redundant calculation", 
	rerender=zz)
yy <- distx(xx1)
princomp(												# reuse computation of distance, only
	xx1, 
	main="a distance computation\ncan be reused too", 
	rerender=yy)
princomp(												# restrict columns analyzed
	xx1, 
	columns= 
		("cow" == columns(xx1, "host_common_name")[,1]))
princomp(
	xx1,												# restrict rows analyzed
	rows= 
		("Carbohydrates" == rows (xx1,"ontology1")[,1]))
princomp(
	xx1, 												# push 3d plot to margins and change persp
	labels="$$project.id",								# with scatterplot3d pars
	map=c(col="host_common_name", pch="samp_store_temp"),
	col=c(Mouse="blue", cow="red", "striped bass"="brown"),
	pch=c("-80"="+",`NA`="x"),
	cex=2,
	angle=20,
	mar=c(1,1,0,0))
princomp(												# label refinement...
	xx1, dim=c(1,2),
	map = c (col="host_common_name", pch="samp_store_temp"),
	col = c (Mouse="brown", cow="red", "striped bass"="blue"),
	pch = c ("-80"="+","NA"="*"),
	cex=2,
	label.font=3, 										# ...italic and
	label.pos=c(1,4,2,2,2,2,4))							# repositioned to stay within box
princomp(
	biom(li4),
	dim=3:1,
	map=c(												# final example with different data
		pch="data.age",
		col="body_site"),
	pch=c(
		"39 ; Year" = 'z',
		"36 ; Year" = 'y',
		"23 ; Year" = 'x'),
	col=c(
		"Teeth surfaces" = "blue"),
	cex=2,
	angle=30,
	box=TRUE,
	box.lty="dashed",
	mar=c(1,1,0,0))

#-----------------------------------------------------------------------------------------
#  image() --- omitted for CRAN
#-----------------------------------------------------------------------------------------
# xx1.log <- transform (xx1, t_Log)
# xx2.log <- transform (xx2, t_Log)
# image(
# 	xx1.log,
# 	margins=c(6,13),
# 	lwid=c(1,1.75), lhei=c(1,10),
# 	cexRow=0.3, cexCol=0.8)
# image(
# 	xx2.log,
# 	margins=c(6,6),
# 	lwid=c(1,2.5), lhei=c(1,10),
# 	cexRow=0.5, cexCol=0.8)
# image(
# 	xx2.log,
# 	margins=c(9,6),
# 	lwid=c(1,2.5), lhei=c(1,10),
# 	cexRow=0.5, cexCol=0.8,
# 	labCol="$$material")
# image(
# 	xx2.log,
# 	margins=c(4,6),
# 	lwid=c(1,2.5), lhei=c(1,10),
# 	cexRow=0.5, cexCol=0.8,
# 	labCol="$$project.id")
#  
# zz <- image (xx1.log)
# image (xx1.log, 											# is this working?
# 	main = "title added without recompute",
# 	margins=c(5,5),
# 	lhei=c(1,3), lwid=c(1,3),
# 	labRow=NA,
# 	rerender=zz)
# 
# image (xx1.log, 											# row subselection
# 	rows = (rows(xx1,"ontology1")[[1]] == "Clustering-based subsystems"),
# 	labRow="$$ontology2",
# 	lwid=c(1,3),
# 	cexRow=0.5,
# 	margins=c(5,10))
# 
# image (xx1.log, columns = c(1,2,4))							# column subselection
# 
# image (xx1.log, labCol=letters[1:7])
# image (xx1.log, labCol = "$$data.age")
# image (xx1.log, labCol=columns(xx1, "data.age") [[1]])		# same as previous
# 
# image (xx1.log, rows=1:20, labRow=1:20)
# image (xx1.log, labRow="$$ontology1")
# image (xx1.log, labRow=rows(xx1, "ontology1")[[1]])			# same as previous
# 
# image(														# no dendrograms
# 	xx2.log,
# 	dendrogram='none',
# 	lwid=c(1,5), lhei=c(1,10),
# 	margins=c(5,7))
