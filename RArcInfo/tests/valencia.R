#1.- Create temporary directory (if needed)
#1.- Extract the E00 file form a ZIP file
#2.- Create an Arc/Info binary coverage
#3.- Create the map

#get current working directory
cwd<-getwd()
#Create tmp directory.
tmpdir<-tempdir()


datadir<-system.file("exampleData",package="RArcInfo")
setwd(datadir)
file.copy(c("valencia.zip", "data_valencia.csv"), tmpdir, overwrite = TRUE)

setwd(tmpdir)

#Convert E00 file to a binary covertage to be imported into R

#Comment this line if the file valencia.e00 already exists
unzip(zipfile="valencia.zip", file="valencia.e00")

#Comments this lines if the binary coverage already exists
library(RArcInfo)
e00toavc("valencia.e00", "valencia")


library(RColorBrewer)
library(RArcInfo)

#Read map data
arcsmuni<-get.arcdata(".", "valencia")
palmuni<-get.paldata(".", "valencia")
bnd.muni<-get.bnddata("info/", "VALENCIA.BND")
patmuni<-get.tabledata("./info", "VALENCIA.PAT")

#Number of polygons
nmuni<-length(palmuni[[1]][[1]])

municipios<-data.frame(1:nmuni, patmuni$"VALENCIA-ID")
names(municipios)<-c("INDEX", "CODMUNICI")


#Datafiles to be used

unemp<-read.table(file="data_valencia.csv", sep=";", 
	dec = ",",skip=1)

unemp<-unemp[,c(1,3)]
names(unemp)<-c("CODMUNICI", "UNEMP")

breaks<-quantile(unemp[,2], c(0,  .025,.2, .8, .975, 1))

unemp<-cbind(unemp, CAT=as.ordered(cut(unemp[,2], breaks, include.lowest = TRUE) ))


#Colors to be used in maps
#colors<-brewer.pal(5, "Oranges")
colors<-brewer.pal(5, "Greens")

ldata<-merge(unemp, municipios, by.x="CODMUNICI", by.y="CODMUNICI")


#Valencia
idx<-(ldata$"CODMUNICI">=46000)
bnd.muni<-c(626679.9, 4250000, 760000, 4460000)


p<-par()
pin<-1.5*p$pin
main<-"Rate of unemployment"

plotpoly(arc=arcsmuni, pal=palmuni, bnd=bnd.muni,
	index=ldata$INDEX[idx], col=colors[ldata$CAT][idx],
	xlab="", ylab="", main=main,
	xaxt="n", yaxt="n", bty="n")

#Set legend
l<-levels(unemp$CAT)
l[1]<-"[0.00,1.26]"
legend(700000, 4460000, fill=colors, 
legend=l, bty="n", cex=1)

setwd(cwd)
