## performance comparison
library(mefa)
library(mefa4)

data(abmibirds)
## NONE: no species was recorded (8)
## SNI: Species Not Identified (216)
## VNA: Variable Not Applicable (6)
## DNC: Did Not Collect (13)
## PNA: Protocol Not Available (0 in this case)

## S3 mefa processing
b3 <- abmibirds
b3 <- b3[!(b3$Scientific.Name %in% c("VNA", "DNC", "PNA")),]
levels(b3$Scientific.Name)[levels(b3$Scientific.Name) %in% c("NONE", "SNI")] <- "zero.pseudo"
b3$Counts <- ifelse(b3$Scientific.Name == "zero.pseudo", 0, 1)
b3$Label <- with(b3, paste(ABMI.Site, Year, Point.Count.Station, sep="_"))
x3 <- b3[!duplicated(b3$Label), c("Label", "ABMI.Site", "Year", "Field.Date", "Point.Count.Station",
    "Wind.Conditions", "Precipitation")]
rownames(x3) <- x3$Label
z3 <- b3[!duplicated(b3$Scientific.Name), c("Common.Name",
    "Scientific.Name", "Taxonomic.Resolution", "Unique.Taxonomic.Identification.Number")]
rownames(z3) <- z3$Scientific.Name
z3 <- z3[z3$Scientific.Name != "zero.pseudo",]
t31 <- system.time(s3 <- suppressWarnings(stcs(b3[,c("Label","Scientific.Name","Counts")])))
t32 <- system.time(m30 <- mefa(s3))
t33 <- system.time(m31 <- mefa(s3, x3, z3))
y30 <- m30$xtab
t34 <- system.time(m32 <- mefa(y30, x3, z3))

## S4 mefa4 processing
b4 <- abmibirds
b4$Label <- with(b4, paste(ABMI.Site, Year, Point.Count.Station, sep="_"))
x4 <- b4[!duplicated(b4$Label), c("Label", "ABMI.Site", "Year", "Field.Date", "Point.Count.Station",
    "Wind.Conditions", "Precipitation")]
rownames(x4) <- x4$Label
z4 <- b4[!duplicated(b4$Scientific.Name), c("Common.Name",
    "Scientific.Name", "Taxonomic.Resolution", "Unique.Taxonomic.Identification.Number")]
rownames(z4) <- z4$Scientific.Name
t41 <- system.time(s4 <- Xtab(~ Label + Scientific.Name, b4, cdrop = c("NONE", "SNI"), 
    subset = !(b4$Scientific.Name %in% c("VNA", "DNC", "PNA")), drop.unused.levels = TRUE))
t42 <- system.time(m40 <- Mefa(s4))
t43 <- system.time(m41 <- Mefa(s4, x4, z4))
y40 <- as.matrix(m40@xtab)
t44 <- system.time(m42 <- Mefa(y40, x4, z4))

res <- cbind("SIZE, *=3"=c("b*"=object.size(b3),
    "s*"=object.size(s3),
    "y*0"=object.size(y30),
    "m*0"=object.size(m30),
    "m*1"=object.size(m31),
    "m*2"=object.size(m32)),
"SIZE, *=4"=c("b*"=object.size(b4),
    "s*"=object.size(s4),
    "y*0"=object.size(y40),
    "m*0"=object.size(m40),
    "m*1"=object.size(m41),
    "m*2"=object.size(m42)),
"TIME, *=3"=c("b*"=NA,
    "s*"=t31[3],
    "y*0"=NA,
    "m*0"=t32[3],
    "m*1"=t33[3],
    "m*2"=t34[3]),
"TIME, *=4"=c("b*"=NA,
    "s*"=t41[3],
    "y*0"=NA,
    "m*0"=t42[3],
    "m*1"=t43[3],
    "m*2"=t44[3]))
(res <- cbind(res, "SIZE"=res[,2]/res[,1], "TIME"=res[,4]/res[,3]))

## debugging

dim(y30)
dim(y40)

setdiff(rownames(y30), rownames(y40))
setdiff(rownames(y40), rownames(y30))

setdiff(colnames(y30), colnames(y40))
setdiff(colnames(y40), colnames(y30))

system.time(xx3 <- aggregate(m31, "ABMI.Site"))
system.time(xx4 <- groupSums(m41, 1, m41@samp$ABMI.Site))
as.Mefa(xx3)
xx4

