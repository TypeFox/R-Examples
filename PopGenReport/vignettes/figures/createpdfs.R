##### File to create included pdfs for the vignette

require(PopGenReport)
data(bilby)

#create the short report on the bilby data set
popgenreport(bilby, mk.counts=T, path.pgr="d:/temp",  fname="PopGenReport_counts")


#create the complete bilby report
popgenreport(bilby, mk.complete=T, path.pgr="d:/temp", fname="PopGenReport_bilby")



#create the map on the tiger occurrences....

tiger.gen <- read.genetable( paste(.libPaths()[1],"/PopGenReport/extdata/tiger.csv",sep="" ), ind=1, pop=2, other.min=3, other.max=6, oneColPerAll=TRUE)


require(rgdal) #load package
xy <- as.matrix(tiger.gen@other$data[,2:1])  #y first then x
#projection from utm 55S to latlong WGS84
latslongs <-project(xy, 
                    "+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs",inv=T) 

#add it to tiger.gen at the right place (again lats)
tiger.gen@other$latlong <- latslongs[,2:1] #again lat and then long
popgenreport(tiger.gen, mk.map=TRUE,mk.counts=FALSE, 
             mapdotcolor="orange", mapdotsize=as.numeric(tiger.gen@other$data$sex), 
             maptype="roadmap", mapdotalpha=0.9,mk.pdf=F, path.pgr="d:/temp", fname="tiger")


#create landgenreport_example.pdf
results1 <-landgenreport(landgen, fric.raster, "D",  path.pgr="d:/temp", fname="landgenreport_example")


#### create platy map
platyfile <-  system.file("extdata/platypus1c.csv",package="PopGenReport")
platy.gen <- read.genetable(platyfile, ind=1, pop=2, lat=3, long=4, other.min=5,
                            other.max=6, oneColPerAll=FALSE, sep="/", ploidy=2)
popgenreport(platy.gen, mk.map=TRUE,mk.counts=FALSE, mapdotcolor="red", 
             maptype="roadmap",mk.pdf=F, path.pgr="d:/temp", fname="platy")