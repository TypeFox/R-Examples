### R code from vignette source 'dicom.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
require("oro.dicom")
require("oro.nifti")
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: DICOM Abdo 01
###################################################
fname <- system.file(file.path("dcm", "Abdo.dcm"), package="oro.dicom")
abdo <- readDICOMFile(fname)
names(abdo)
head(abdo$hdr)
tail(abdo$hdr)


###################################################
### code chunk number 3: DICOM Abdo 02
###################################################
extractHeader(abdo$hdr, "BitsAllocated")
extractHeader(abdo$hdr, "Rows")
extractHeader(abdo$hdr, "Columns")


###################################################
### code chunk number 4: DICOM HK40 00
###################################################
load(system.file("hk-40/hk40.RData", package="oro.dicom"))


###################################################
### code chunk number 5: DICOM HK40 01 (eval = FALSE)
###################################################
## fname <- system.file("hk-40", package="oro.dicom")
## data(dicom.dic)
## hk40 <- readDICOM(fname)


###################################################
### code chunk number 6: DICOM HK40 01.1
###################################################
unlist(lapply(hk40, length))


###################################################
### code chunk number 7: DICOM HK40 02
###################################################
hk40.info <- dicomTable(hk40$hdr)
write.csv(hk40.info, file="hk40_header.csv")
sliceloc.col <- which(hk40$hdr[[1]]$name == "SliceLocation") 
sliceLocation <- as.numeric(hk40.info[, sliceloc.col])
head(sliceLocation)
head(diff(sliceLocation))
unique(extractHeader(hk40$hdr, "SliceThickness"))


###################################################
### code chunk number 8: DICOM HK40 03
###################################################
head(extractHeader(hk40$hdr, "SliceLocation"))
modality <- extractHeader(hk40$hdr, "Modality", numeric=FALSE)
head(matchHeader(modality, "mr"))
(seriesTime <- extractHeader(hk40$hdr, "SeriesTime", numeric=FALSE))
str2time(seriesTime)


###################################################
### code chunk number 9: abdo-png
###################################################
jpeg(filename="dicom_abdo.jpeg", width=480, height=480, quality=95, bg="black")
par(mar=rep(0,4))


###################################################
### code chunk number 10: abdo-image
###################################################
image(t(abdo$img), col=grey(0:64/64), axes=FALSE, xlab="", ylab="")


###################################################
### code chunk number 11: abdo-dev.off
###################################################
dev.off()


###################################################
### code chunk number 12: DICOM Abdo 03
###################################################
extractHeader(abdo$hdr, "Manufacturer", numeric=FALSE)
extractHeader(abdo$hdr, "RepetitionTime")
extractHeader(abdo$hdr, "EchoTime")


###################################################
### code chunk number 13: DICOM Siemens 01
###################################################
fname <- system.file(file.path("dcm", "MR-sonata-3D-as-Tile.dcm"),
                     package="oro.dicom")
dcm <- readDICOMFile(fname)
dim(dcm$img)
dcmImage <- create3D(dcm, mosaic=TRUE)
dim(dcmImage)


###################################################
### code chunk number 14: DICOM Siemens 02
###################################################
dcmNifti <- dicom2nifti(dcm, mosaic=TRUE)
jpeg(filename="dcmImage.jpeg", width=2*480, height=2*480, quality=95, bg="black")
image(t(dcm$img), col=grey(0:64/64), zlim=c(16, 1024), axes=FALSE, 
      xlab="", ylab="")
dev.off()
jpeg(filename="dcmNifti.jpeg", width=480, height=480, quality=95, bg="black")
image(dcmNifti, zlim=c(16, 1024))
dev.off()


###################################################
### code chunk number 15: DICOM2NIFTI HK40 01
###################################################
dput(formals(dicom2nifti))
(hk40n <- dicom2nifti(hk40))


###################################################
### code chunk number 16: DICOM2NIFTI HK40 02 (eval = FALSE)
###################################################
## image(hk40n)
## orthographic(hk40n, col.crosshairs="green")


###################################################
### code chunk number 17: DICOM2NIFTI HK40 03
###################################################
jpeg("hk40n_image.jpeg", width=480, height=480, quality=95, bg="black")
image(hk40n, zlim=c(0,1024))
dev.off()
jpeg("hk40n_orthographic.jpeg", width=480, height=480, quality=95, bg="black")
orthographic(hk40n, zlim=c(0,1024), col.crosshairs="green")
dev.off()


###################################################
### code chunk number 18: DICOM2NIFTI HK40 04
###################################################
(hk40n <- dicom2nifti(hk40, DIM=4))


###################################################
### code chunk number 19: RIDER Neuro MRI (eval = FALSE)
###################################################
## subject <- "1086100996"
## DCM <- readDICOM(subject, verbose=TRUE)
## seriesInstanceUID <- extractHeader(DCM$hdr, "SeriesInstanceUID", FALSE)
## for (uid in unique(seriesInstanceUID)) {
##   index <- which(unlist(lapply(DCM$hdr, function(x) uid %in% x$value)))
##   uid.dcm <- list(hdr=DCM$hdr[index], img=DCM$img[index])
##   patientsName <- extractHeader(uid.dcm$hdr, "PatientsName", FALSE)
##   studyDate <- extractHeader(uid.dcm$hdr, "StudyDate", FALSE)
##   seriesDescription <- extractHeader(uid.dcm$hdr, "SeriesDescription", FALSE)
##   fname <- paste(gsub("[^0-9A-Za-z]", "", 
##                       unique(c(patientsName, studyDate, seriesDescription))), 
##                  collapse="_")
##   cat("##  ", fname, fill=TRUE)
##   if (gsub("[^0-9A-Za-z]", "", unique(seriesDescription)) == "axtensor") {
##     D <- 4
##     reslice <- FALSE
##   } else {
##     D <- 3
##     reslice <- TRUE
##   }
##   uid.nifti <- dicom2nifti(uid.dcm, DIM=D, reslice=reslice,
##                            descrip=c("PatientID", "SeriesDescription"))
##   writeNIfTI(uid.nifti, fname)
## }


