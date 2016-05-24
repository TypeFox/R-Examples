##    DATforDCEMRI: a Deconvolution Analysis Tool for DCE-MRI
##    Copyright 2013 Genentech, Inc.
##
##    This package is distributed under the
##    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License
##    at http://creativecommons.org/licenses/by-nc-sa/3.0/us/
##
##    For questions or comments, please contact
##    Gregory Z. Ferl, Ph.D.
##    Genentech, Inc.
##    Development Sciences
##    1 DNA Way, Mail stop 463A
##    South San Francisco, CA, United States of America
##    E-mail: ferl.gregory@gene.com



DAT <-
function(file="nodata", slice=0, vp=0, border=20, maxCt=0.66, parameter.plot="AUCMRT", cutoff.map=0.85, range.map=1.5, export.matlab=FALSE, batch.mode=FALSE, alpha.AIF=c(0,0,2000), correct.trunc=TRUE, kep.nom=0.5, ...){

    ## Workaround for "DAT: no visible binding for global variable 'DAT.simData'" NOTE
    ## Lazy Load instead
    ## DAT.simData <- NULL

if(export.matlab==TRUE)
  require(R.matlab)
require(grid)
require(graphics)
require(locfit)
require(matlab)
require(xtable)
require(akima)
require(R.oo)
require(R.methodsS3)
require(lattice)

AUMC <- function(AUMC.median, h.median, irf_time_vec, r){
    #AUMC.median <- AUMC.median + 0.5*(h.median[r] + h.median[r+1])*0.5*(irf_time_vec[r] + irf_time_vec[r+1])
     AUMC.median <- AUMC.median + 0.5*(h.median[r]*irf_time_vec[r] + h.median[r+1]*irf_time_vec[r+1])
}

#print(AUMC)

DATrun <- function(file, slice, vp, border, maxCt, parameter.plot, cutoff.map, range.map, export.matlab, batch.mode, alpha.AIF, correct.trunc){

DAT.version <- 0.55
ptm_total <- proc.time()[3]


###############################################################
####BASIC DECONVOLUTION FUNCTION###############################
###############################################################
calch <- function(u, y, TIME_trunc){
###########TRY ALTERNATE PREPLOT PARAMETER HERE#################
locfit_y <- preplot(y, newdata=0:max(TIME_trunc))
#locfit_y <- preplot(y, newdata=TIME_trunc)
y_smooth <- locfit_y$fit

n<-length(u)
A<-matrix(0,nrow=n,ncol=n)
ind<-row(A)-col(A)
ind[ind<0] <- (-1)
ind<-ind+2
A <- matrix(c(0,u)[ind],nrow=n,ncol=n)
h <- solve(A,y_smooth)
calch_out <- list(h, y_smooth)
names(calch_out) <-c("h", "y_smooth")
return(calch_out)
}


###################################################
###PACKAGE INFORMATION#############################
###################################################

cat("\n")
cat("#######################################################################", "\n")
cat("#####-- DATforDCEMRI: a Deconvolution Analysis Tool for DCE-MRI --#####", "\n")
cat("#####------------------- R package version",DAT.version,"------------------#####", "\n")
cat("#####--------------- Copyright 2013 Genentech, Inc. --------------#####", "\n")
cat("#######################################################################", "\n")
cat("#####--------- For questions and comments please contact ---------#####", "\n")
cat("#####---------------------- Gregory Z. Ferl ----------------------#####", "\n")
cat("#####------------------- ferl.gregory@gene.com -------------------#####", "\n")
cat("#######################################################################", "\n")
cat("#####------ DATforDCEMRI comes with ABSOLUTELY NO WARRANTY. ------#####", "\n")
cat("#####------- This is free software, and you are welcome to -------#####", "\n")
cat("#####--------- redistribute it under certain conditions. ---------#####", "\n")
cat("#####-------------------- See Creative Commons -------------------#####", "\n")
cat("#####------ Attribution-NonCommercial-ShareAlike 3.0 License -----#####", "\n")
cat("#####-- at http://creativecommons.org/licenses/by-nc-sa/3.0/us/ --#####", "\n")
cat("#####------------------------ for details. -----------------------#####", "\n")
cat("#######################################################################", "\n")
cat("\n")


###################################################
###LOAD THE ORIGINAL MATLAB FILE INTO R############
###################################################

###IMPORT MATLAB FILE INTO R#######################
filea <- strsplit(file, split="/")[[1]]
fileb <- filea[length(filea)]

cat("loading", fileb, "into R...", "\n")
ptm <- proc.time()[3]

load(file)
data <- dcemri.data
rm(dcemri.data)
cat("done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")

if(length(unique(names(data)=="status")) == 2)
  status.maps.included <- TRUE
if(length(unique(names(data)=="status")) == 1)
  status.maps.included <- FALSE

if(length(names(data)) < 8)
  mat.original <- TRUE
if(length(names(data)) >= 8)
  mat.original <- FALSE

if(mat.original == TRUE){

if(length(dim(data$maskROI))==3){
  for(i in length(data$maskROI[1,1,]):1){
    if(max(data$maskROI[,,i])==1)
      roi_start <- i
  }
  for(i in 1:length(data$maskROI[1,1,])){
    if(max(data$maskROI[,,i])==1)
      roi_end <- i
  }
  if(slice=="" || slice < roi_start || slice > roi_end)
    stop(paste("Please specify a slice number between ", roi_start, " and ", roi_end, ".", sep=""))
}



###EXTRACT THE TIME VECTOR AND CONVERT TO MINUTES########################
vector.times <- as.vector(data$vectorTimes/60)

###EXTRACT CONTRAST AGENT CONCENTRATIONS FOR SLICE OF INTEREST###########
map_cc <- data$mapCC

if(length(dim(data$mapCC))==4){
  map_cc_slice <- map_cc[,,slice,]
  cat("extracting and processing slice", slice, "for analysis...", "\n")
  ptm <- proc.time()[3]
}

if(length(dim(data$mapCC))==3){
  map_cc_slice <- map_cc
  cat("processing slice for analysis...", "\n")
  ptm <- proc.time()[3]
}

if(length(vector.times)!=length(map_cc_slice[1,1,]))
  stop("The length of the time vector does not match the length of the contrast agent concentration vector")

nt <- length(map_cc_slice[1,1,])
ny <- length(map_cc_slice[,1,1])
nx <- length(map_cc_slice[1,,1])

for(i in 1:nt)
 map_cc_slice[,,i] <- rot90(map_cc_slice[,,i],3)


###EXTRACT ROI MASK FOR SLICE OF INTEREST################################
if(length(dim(data$maskROI))==3)
  mask.roi <- data$maskROI[,,slice]

if(length(dim(data$maskROI))==2)
  mask.roi <- data$maskROI

mask.roi <- rot90(mask.roi,3)

if(is.finite(median(mask.roi))==FALSE)
  stop("Your ROI mask contains nonnumeric elements; voxels within the ROI should have a value of ``1'' and all other voxels should have a value of ``0''.")

if(length(mask.roi)==length(mask.roi[mask.roi==0]))
  stop("Your ROI mask is composed entirely of zeroes; voxels within the ROI should have a value of ``1'' and all other voxels should have a value of ``0''.")

if(length(mask.roi)==length(mask.roi[mask.roi==1]))
  warning("Your ROI mask is composed entirely of ones.")

if(max(mask.roi)>1 || min(mask.roi)<0)
  stop("Your ROI mask contains values greater than one and/or less than zero; voxels within the ROI should have a value of ``1'' and all other voxels should have a value of ``0''.")


###EXTRACT THE AIF VECTOR FROM THE DATA FILE#########################
AIF <- as.vector(data$vectorAIF)

rm(data)
cat("done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")

}



##############################################################
####GENERATE A UNIQUE FILENAME BASED ON THE CURRENT DATE######
##############################################################
filename <- strsplit(file, split="/")[[1]]

filename <- filename[length(filename)]
filename <- strsplit(filename, split=".mat")[[1]]

if(slice==0)
  filename2 <- ""
if(slice!=0)
  filename2 <- paste("s", slice, sep="")

DATE <- date()
if(length(strsplit(DATE,split="  ")[[1]])==2){
  full_date <- strsplit(DATE,split="  ")[[1]]
  full_date_1 <- full_date[1]
  full_date_1 <- strsplit(full_date_1, split=" ")[[1]]
  full_date_2 <- full_date[2]
  full_date_2 <- strsplit(full_date_2, split=" ")[[1]]
  month <- full_date_1[2]
  day   <- full_date_2[1]
  time  <- full_date_2[2]
  year  <- full_date_2[3]
  }
if(length(strsplit(DATE,split="  ")[[1]])==1){
  full_date <- strsplit(DATE,split=" ")[[1]]
  month <- full_date[2]
  day <- full_date[3]
  time <- full_date[4]
  year <- full_date[5]
  }
year <- strsplit(year, split="")[[1]]
year <- paste(year[3],year[4],sep="")
time_concat <- strsplit(time, split=":")[[1]]
time_concat <- paste(paste(time_concat[1], time_concat[2], sep=""), time_concat[3], sep="")
DATE <- paste(paste(paste(paste(day,month,sep=""), year, sep=""), "-", sep=""), time_concat, sep="")

filename3 <- paste("DAT_", filename, "_", filename2, "_", DATE, ".mat", sep="")



if(mat.original == TRUE){
cat("performing deconvolution analysis on entire FOV...", "\n")
ptm <- proc.time()[3]
artery_data <- data.frame(vector.times*60, AIF)
names(artery_data) <- c("TIME", "ARTERY")
data_artery_peak <- subset(artery_data, artery_data$ARTERY == max(artery_data$ARTERY))
data_remove_artery_prepeak <- subset(artery_data, artery_data$TIME >= data_artery_peak$TIME)
frames_to_peak <- length(artery_data[,1]) - length(data_remove_artery_prepeak[,1]) + 1
TIME   <- data_remove_artery_prepeak$TIME
ARTERY <- data_remove_artery_prepeak$ARTERY
TIME_trunc <- TIME[seq(1,length(TIME)-1,by=1)]
TIME_trunc <- TIME_trunc - TIME_trunc[1]
ARTERY_trunc <- ARTERY[seq(1,length(ARTERY)-1,by=1)]

ARTERY_smooth <- locfit.robust(ARTERY_trunc~TIME_trunc, acri="cp", alpha=alpha.AIF)
AIF_smooth <- ARTERY_smooth

###########TRY ALTERNATE PREPLOT PARAMETER HERE#################
#locfit_u <- preplot(AIF_smooth, newdata=TIME_trunc)
locfit_u <- preplot(AIF_smooth, newdata=0:max(TIME_trunc))
u_smooth <- locfit_u$fit

Tmax <- max(TIME_trunc)

#####################################################################################
###########Perform deconvolution analysis on median voxel-wise curves################
#####################################################################################
TUMOR.median <- 1:length(map_cc_slice[1,1,])
for(i in 1:length(map_cc_slice[1,1,]))
  TUMOR.median[i] <- median(map_cc_slice[,,i][mask.roi==1], na.rm=TRUE)

TUMOR.median <- TUMOR.median[seq(frames_to_peak,length(TUMOR.median),by=1)]

TUMOR.median_corr <- TUMOR.median - vp*ARTERY
TUMOR.median_corr_shifted <- TUMOR.median_corr[seq(2,length(TUMOR.median_corr),by=1)]
TUMOR.median_smooth <- locfit.robust(TUMOR.median_corr_shifted~TIME_trunc, acri="cp")

calch.out <- calch(u_smooth, TUMOR.median_smooth, TIME_trunc)
h.median <- calch.out$h
tumor_smooth_median <- calch.out$y_smooth

#n <- length(h.median)
#n2 <- Tmax
#interval <- n2/(n-1)
#irf_time_vec <- seq(0, n2, by=interval)

irf_time_vec <- 0:max(TIME_trunc)
n <- length(h.median)

AUC.median <- 0
AUMC.median <- 0
  for(r in 1:(n-1)){
    h_sum <- h.median[r] + h.median[r+1]
    t_sum <- irf_time_vec[r] + irf_time_vec[r+1]
    AUC.median <- AUC.median + 0.5*h_sum
    AUMC.median <- AUMC.median + 0.5*(h.median[r]*irf_time_vec[r] + h.median[r+1]*irf_time_vec[r+1])
  }
AUCMRT.median <- AUC.median/(AUMC.median/AUC.median)*60

#####################################################################################
#####################################################################################




assign("map.AUC", array(rep(0,nx*ny), dim=c(nx,ny)))
assign("map.AUCMRT", array(rep(0,nx*ny), dim=c(nx,ny)))
assign("tumor_data_array", array(rep(0,nx*ny*length(ARTERY_trunc)), dim=c(nx,ny,length(ARTERY_trunc))))
assign("tumor_smooth_array", array(rep(0,nx*ny*length(u_smooth)), dim=c(nx,ny,length(u_smooth))))
assign("map.h", array(rep(0,nx*ny*length(u_smooth)), dim=c(nx,ny,length(u_smooth))))

for(i in border:(nx-border)){
  for(j in border:(ny-border)){
      TUMOR <- map_cc_slice[i,j,]
      TUMOR <- TUMOR[seq(frames_to_peak,length(TUMOR),by=1)]
        if(is.finite(mean(TUMOR))!= FALSE){
          if(max(TUMOR) < maxCt*max(ARTERY) & max(TUMOR)>0){
            TUMOR_corr <- TUMOR - vp*ARTERY
            TUMOR_corr_shifted <- TUMOR_corr[seq(2,length(TUMOR_corr),by=1)]
            tumor_data_array[i,j,] <- TUMOR_corr_shifted

            TUMOR_smooth <- locfit.robust(TUMOR_corr_shifted~TIME_trunc, acri="cp")

            calch.out <- calch(u_smooth, TUMOR_smooth, TIME_trunc)
            h <- calch.out$h
            tumor_smooth_array[i,j,] <- calch.out$y_smooth

            #if(abs(min(h)) < 0.1*abs(max(h[2:length(h)])) || min(h) >= 0){
              ####CALCULATE AUC AND MRT OF A CONC_TIME CURVE#################
              #n <- length(h)
              #n2 <- Tmax
              #interval <- n2/(n-1)
              #irf_time_vec <- seq(0, n2, by=interval)
              #h_wtime <- cbind(irf_time_vec, h)

              n <- length(h)
              irf_time_vec <- 0:max(TIME_trunc)

              AUC <- 0
              AUMC <- 0

              ###TEST THIS TIME VEC#######
              #irf_time_vec <- 0:max(TIME_trunc)

              for(r in 1:(n-1)){
                #h_sum <- h[r] + h[r+1]
                #t_sum <- irf_time_vec[r] + irf_time_vec[r+1]
                #AUC <- AUC + 0.5*h_sum
                #AUMC <- AUMC + 0.25*h_sum*t_sum
                AUC <- AUC + 0.5*(h[r] + h[r+1])
                AUMC <- AUMC + 0.25*(h[r] + h[r+1])*(irf_time_vec[r] + irf_time_vec[r+1])
              }

             map.AUC[i,j] <- AUC
             map.AUCMRT[i,j] <- AUC/(AUMC/AUC)*60
             map.h[i,j,] <- h

            #}
          }
        }



  }

}


######CORRECT FOR TRUNCATION ERRORS - ROUGH APPROXIMATE ONLY#################
if(correct.trunc==TRUE){
  #kep_nom <- median(map.AUCMRT[map.AUCMRT>0])/median(map.AUC[map.AUC>0])
  t_scan <- max(TIME_trunc)/60
  ve_trunc_error <- 1-exp(-kep.nom*t_scan)
  Ktrans_trunc_error <- (1-exp(-kep.nom*t_scan))^2/(1-(1+kep.nom*t_scan)*exp(-kep.nom*t_scan))
  map.AUC <- map.AUC / ve_trunc_error
  map.AUCMRT <- map.AUCMRT / Ktrans_trunc_error

  AUC.median <- AUC.median / ve_trunc_error
  AUMC.median <- AUMC.median / Ktrans_trunc_error
}

cat("..done in", format((proc.time()[3]-ptm)/60, digits=2), "minutes.", "\n")
}






if(mat.original==FALSE){


##############################################################
####LOAD MATLAB FILE PREVIOUSLY GENERATED BY THIS SCRIPT######
##############################################################
if(mat.original == FALSE){
load(file)
data <- dcemri.data
rm(dcemri.data)
args <- data$args
plotParams <- data$plotParams
nx <- plotParams$nx
ny <- plotParams$ny
nt <- plotParams$nt
file <- args$file
slice <- args$slice
vp <- args$vp
border <- args$border
maxCt <- args$maxCt
#parameter.plot <- args$parameterplot
range.map <- args$rangemap
cutoff.map <- args$cutoffmap
export.matlab <- args$exportmatlab
batch.mode <- args$batchmode
alpha.AIF <- args$alphaAIF
correct.trunc <- args$correcttrunc
kep.nom <- args$kepnom
map_cc_slice <- data$mapCC
vector.times <- data$vectorTimes
AIF <- data$vectorAIF
mask.roi <- data$maskROI
map.AUC <- data$mapAUC
map.AUCMRT <- data$mapAUCMRT
map.h <- data$mapIRF
tumor_data_array <- data$mapCCtransformed
tumor_smooth_array <- data$mapCCsmoothed
AIF_smooth<- data$vectorAIFsmoothed
ARTERY_trunc <- data$vectorAIFtrunc
TIME_trunc<- data$vectorTimesTrunc
tumor_smooth_median <- data$CCmedianSmoothed
h.median <- data$IRFmedian
AUC.median <- data$AUCmedian
AUCMRT.median <- data$AUCMRTmedian
TUMOR.median_corr_shifted <- data$CCmedianTransformed
filename3 <- data$filenameS
rm(data)

h.median <- h.median[2:length(h.median)]
map.h <- map.h[,,2:length(h.median)]
}
######################################################################

if(parameter.plot=="AUCMRT"){
  median.value <- format(median(map.AUCMRT[mask.roi==1], na.rm=TRUE), digits=3)

  map_ul <- map.AUCMRT[map.AUCMRT>0]
  map_ul <- sort(map_ul)
  map_ul <- range.map*(max(map_ul[1:length(map_ul)*cutoff.map]))
  map.AUCMRT[map.AUCMRT<0] <- 0
  map.AUCMRT[map.AUCMRT>=map_ul*0.99] <- map_ul*0.99

  map.deconv.toplot <- map.AUCMRT

  plot.title <- "AUC divided by the Mean Residence Time of the Impulse Response Function"
  plot.units <- "1/min"
}

if(parameter.plot=="AUC"){
  median.value <- format(median(map.AUC[mask.roi==1], na.rm=TRUE), digits=3)


  map_ul <- map.AUC[map.AUC>0]
  map_ul <- sort(map_ul)
  map_ul <- range.map*(max(map_ul[1:length(map_ul)*cutoff.map]))
  map.AUC[map.AUC<0] <- 0
  map.AUC[map.AUC>=map_ul*0.99] <- map_ul*0.99

  map.deconv.toplot <- map.AUC

  plot.title <- "AUC of the Impulse Response Function"
  plot.units <- "dimensionless"

}

############################################
###DRAW THE ROI ON THE DECONVOLUTION MAP####
############################################
for(x in 1:nx){
  for(y in 1:ny){
    if(mask.roi[x,y]==1){
      map.deconv.toplot[x,y] <- NA
      y <-1
      break
    }
  }
y <- 1
}

for(x in nx:1){
  for(y in ny:1){
    if(mask.roi[x,y]==1){
      map.deconv.toplot[x,y] <- NA
      y <- ny
      break
    }
  }
y <- ny
}

for(y in 1:ny){
  for(x in 1:nx){
    if(mask.roi[x,y]==1){
      map.deconv.toplot[x,y] <- NA
      x <-1
      break
    }
  }
x <- 1
}

for(y in ny:1){
  for(x in nx:1){
    if(mask.roi[x,y]==1){
      map.deconv.toplot[x,y] <- NA
      x <- nx
      break
    }
  }
x <- nx
}

dev.new(width=7.675, height=7.675, xpos=240, ypos=0)
image(map.deconv.toplot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), main=paste(plot.title), zlim=c(0,map_ul), cex.main=0.75)
image(map.deconv.toplot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), main=paste(plot.title), zlim=c(0,map_ul), cex.main=0.75)
legend_deconv <- seq(0, map_ul, by=0.001)
dim(legend_deconv) <- c(1,length(legend_deconv))

dev.new(width=2.5, height=7.675, xpos=0, ypos=0)
image(legend_deconv, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xaxt="n", yaxt="n", xlab="0", ylab=plot.units, main=format(map_ul, digits=3), cex.lab=1.25, cex.main=1.05)




}

###########################################################################################################################################
########NOTE THAT CALCULATION OF x_max, x_min, y_max, y_min IS PERFORMED WHEN mat.original=TRUE AND mat.original=FALSE#####################
###########################################################################################################################################



####ZOOM IN ON TUMOR ROI#####
x_min <- 0
x_max <- nx
y_min <- 0
y_max <- ny

    for(xx in 1:nx){
      if(sum(mask.roi[xx,], na.rm=TRUE)>0){
        x_min <- xx - 3
        break
      }
      if(sum(mask.roi[xx,], na.rm=TRUE)==0)
        xx <- xx+1
     }

    for(y in 1:ny){
      if(sum(mask.roi[,y], na.rm=TRUE)>0){
        y_min <- y - 3
        break
      }
      if(sum(mask.roi[,y], na.rm=TRUE)==0)
        y <- y+1
     }

    for(xx in nx:0){
      if(sum(mask.roi[xx,], na.rm=TRUE)>0){
        x_max <- xx + 3
        break
      }
      if(sum(mask.roi[xx,], na.rm=TRUE)==0)
        xx <- xx-1
     }

    for(y in ny:0){
      if(sum(mask.roi[,y], na.rm=TRUE)>0){
        y_max <- y + 3
        break
      }
      if(sum(mask.roi[,y], na.rm=TRUE)==0)
        y <- y-1
     }


if(parameter.plot == "AUC")
  MAP <- map.AUC
if(parameter.plot == "AUCMRT")
  MAP <- map.AUCMRT

MAP_ul_s1 <- MAP[MAP>0]
MAP_ul_s1 <- sort(MAP_ul_s1)
MAP_ul <- range.map*(max(MAP_ul_s1[1:length(MAP_ul_s1)*cutoff.map]))


#if(mat.original==FALSE || batch.mode==FALSE){
if(mat.original==FALSE){

if(parameter.plot == "AUC")
  MAP <- map.AUC
if(parameter.plot == "AUCMRT")
  MAP <- map.AUCMRT

if(parameter.plot == "AUC")
  param.median <- median(map.AUC[map.AUC>0], na.rm=TRUE)

if(parameter.plot == "AUCMRT")
  param.median <- median(map.AUCMRT[map.AUCMRT>0], na.rm=TRUE)



####CALCULATE RANGE OF COLOR BAR (LEGEND)#####

MAP_for_plot <- MAP

####REPLACE ANY NEGATIVE NUMBERS WITH ZEROS######
MAP_for_plot[MAP_for_plot<0] <- 0

####TRUNCATE HIGH OUTLIER VALUES#####
MAP_for_plot[MAP_for_plot>=MAP_ul*0.99] <- MAP_ul*0.99

dev.new(width=7.675, height=7.675, xpos=964, ypos=0)

image(MAP_for_plot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), zlim=c(0,MAP_ul), main=paste(filename, "       slice ", args$slice, "       median ", parameter.plot, " (within ROI) = " , median.value, sep=""), cex.main=0.75)
image(MAP_for_plot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), zlim=c(0,MAP_ul), main=paste(filename, "       slice ", args$slice, "       median ", parameter.plot, " (within ROI) = " , median.value, sep=""), cex.main=0.75)

text(x_min/nx+(x_max/nx-x_min/nx)*0.96, y_min/ny+(y_max/ny-y_min/ny)*0.98, "close", col="green")
text(x_min/nx+(x_max/nx-x_min/nx)*0.08, y_min/ny+(y_max/ny-y_min/ny)*0.98, "print to PDF", col="green")
text(x_min/nx+(x_max/nx-x_min/nx)*0.135, y_min/ny+(y_max/ny-y_min/ny)*0.02, "display median IRF", col="green")
text(x_min/nx+(x_max/nx-x_min/nx)*0.73, y_min/ny+(y_max/ny-y_min/ny)*0.02, paste("DAT for DCE-MRI:", DAT.version), col="green")
legend <- seq(0,MAP_ul,by=0.001)
dim(legend) <- c(1,length(legend))



#########VOXEL LOCATOR####################
dev.set(4)

inf <- 1
newplot <- 1

legend_count<-1
legend_labels <- 1:1000
legend_matrix <- matrix(0, ncol=nx, nrow=ny)

cat("---", "\n")

while(inf == 1){

z <- locator(1, type="o", col="green")

xx <- round(z$x*(nx-1)+1)
yy <- round(z$y*(ny-1)+1)

if(legend_matrix[xx,yy]==0){
legend(z$x, z$y, legend_labels[legend_count], col="green", text.col="green")
legend_matrix[xx,yy] <- 1
cat("point", legend_count, "\n")
cat("x =", xx, "\n")
cat("y =", yy, "\n")
cat("AUC =", format(map.AUC[xx,yy],digits=3), "; ")
cat("AUC/MRT =", format(map.AUCMRT[xx,yy], digits=3), "\n")
cat("---", "\n")
legend_count <- legend_count+1
}

xdim <- x_max-x_min
ydim <- y_max-y_min

if(xx > (x_max-0.1*xdim) && yy > (y_max-0.03*ydim)){
  graphics.off()
  cat("---", "\n")
  cat("session ended", "\n")
  cat("---", "\n")
  break()
}

if(xx > (x_max-0.5*xdim) && yy < (y_min+0.06*ydim)){
cat("---", "\n")
cat("    ###########################################################################", "\n")
cat("    DATforDCEMRI: a Deconvolution Analysis Tool for DCE-MRI", "\n")
cat("    Copyright 2011 Genentech, Inc.", "\n")
cat("\n")
cat("    This program is free software; you can redistribute it and/or modify", "\n")
cat("    it under the terms of the GNU General Public License as published by", "\n")
cat("    the Free Software Foundation; either version 2 of the License, or", "\n")
cat("    (at your option) any later version.", "\n")
cat("\n")
cat("    This program is distributed in the hope that it will be useful,", "\n")
cat("    but WITHOUT ANY WARRANTY; without even the implied warranty of", "\n")
cat("    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the", "\n")
cat("    GNU General Public License for more details.", "\n")
cat("\n")
cat("    You should have received a copy of the GNU General Public License", "\n")
cat("    along with this program; if not, write to the Free Software", "\n")
cat("    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA", "\n")
cat("\n")
cat("    For questions or comments, please contact", "\n")
cat("    Gregory Z. Ferl, Ph.D.", "\n")
cat("    Genentech, Inc.", "\n")
cat("    Development Sciences", "\n")
cat("    1 DNA Way, Mail stop 463A", "\n")
cat("    South San Francisco, CA, United States of America", "\n")
cat("    E-mail: ferl.gregory@gene.com", "\n")
cat("    ###########################################################################", "\n")
cat("---", "\n")
}

xx_old <- xx
yy_old <- yy




######PLOT DATA AND SIMULATION############

conc <- 1:nt

for(i in 1:nt)
  conc[i] <- map_cc_slice[xx,yy,i]



if(xx < (x_min + 0.12*xdim) && yy > (y_max-0.03*ydim)){

############PRINT SUMMARY FIGURES TO PDF FILE################################################################
pdf(file=paste(strsplit(filename3,split=".mat")[[1]], ".pdf", sep=""), height=12, width=18)
#pdf(file=sub(".RData", ".pdf", file), height=12, width=18)

layout(matrix(c(1,2,3,4,5,6,7,8), 2, 4, byrow=TRUE), widths=c(1.5,5.5,5.5,5.5))
par(omi=c(0.15, 0.15, 0.15, 0.15))

image(legend_deconv, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xaxt="n", yaxt="n", xlab="0", ylab=plot.units, main=format(map_ul, digits=3), cex.lab=1.5, cex.main=1.25, cex.axis=1.5)

image(map.deconv.toplot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), main=plot.title, zlim=c(0,map_ul), cex.main=1.25, cex.axis=1.5, cex.lab=1.5)

image(MAP_for_plot, col=palette(colorRampPalette(c("black", "red","yellow"))(1000)), xlim=c(x_min/nx,x_max/nx), ylim=c(y_min/ny,y_max/ny), zlim=c(0,MAP_ul), cex.main=1.25, cex.axis=1.25, main=paste("median ", parameter.plot, " (within ROI) = " , median.value, sep=""))

plot(0.5, 0.5, xaxt="n", yaxt="n", ylim=c(0,1), xlim=c(0,1), xlab="", ylab="", col="white", main="Function arguments for DAT()", cex.main=1.5)
text(0.5, 0.97,  DAT.version, cex=2.00)
text(0.5, 0.88,  paste("file =", file), cex=2.00)
text(0.5, 0.81,  paste("slice =", slice), cex=2.00)
text(0.5, 0.74,  paste("vp =", vp), cex=2.00)
text(0.5, 0.67,  paste("border =", border), cex=2.00)
text(0.5, 0.60,  paste("maxCt =", maxCt), cex=2.00)
text(0.5, 0.53,  paste("parameter.plot =", parameter.plot), cex=2.00)
text(0.5, 0.46,  paste("range.map =", range.map), cex=2.00)
text(0.5, 0.39,  paste("cutoff.map =", cutoff.map), cex=2.00)
text(0.5, 0.32,  paste("export.matlab =", export.matlab), cex=2.00)
text(0.5, 0.25,  paste("batch.mode =", batch.mode), cex=2.00)
text(0.5, 0.18,  paste("alpha.AIF = c(", alpha.AIF[1],", ", alpha.AIF[2],", ", alpha.AIF[3], ")", sep=""), cex=2.00)
text(0.5, 0.11,  paste("correct.trunc =", correct.trunc), cex=2.00)


frame()
plot(TIME_trunc, ARTERY_trunc, xlab="min", ylim=c(0, max(ARTERY_trunc)), ylab="micro-mol/L", main="Arterial Input Function with smoothed curve", cex=3, cex.lab=1.5, cex.main=1.5, cex.axis=1.5, lwd=2)
lines(AIF_smooth, lwd=3)

plot(TIME_trunc, TUMOR.median_corr_shifted, xlab="min", ylim=c(0, max(TUMOR.median_corr_shifted)), ylab="micro-mol/L", main="Median contrast agent data with smoothed curve", cex=3, cex.lab=1.5, cex.main=1.5, cex.axis=1.5, lwd=2)
lines(tumor_smooth_median, lwd=3)

plot(h.median, xlab="sec", ylim=c(0, max(h.median)), ylab=expression(sec^-1), main="Median Impulse Response Function", col="white", cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
lines(h.median, lwd=3, col="goldenrod")
text(max(TIME_trunc)/2, max(h.median)*0.93, paste(paste("AUC =", format(AUC.median, digits=3))), cex=1.75)
text(max(TIME_trunc)/2, max(h.median)*0.86, paste("AUC/MRT =", format(AUCMRT.median, digits=3), "1/min"), cex=1.75)

dev.off()

cat("image printed to ", sub(".RData", ".pdf", paste(strsplit(filename3,split=".mat")[[1]], ".pdf", sep="")), ".", "\n", sep="")
cat("---", "\n")
#############################################################################################################

}

############PLOT DECONVOLUTION CURVES#################

if(newplot==1){
dev.new(width=5.95, height=3.95, xpos=0, ypos=734)
plot(TIME_trunc, ARTERY_trunc, xlab="min", ylim=c(0, max(ARTERY_trunc)), ylab="micro-mol/L", main="AIF (nonparametric smoother)")
lines(AIF_smooth)
dev.new(width=5.95, height=3.95, xpos=564, ypos=734)
dev.new(width=5.95, height=3.95, xpos=1124, ypos=734)
newplot <- 2
}



dev.set(6)
plot(TIME_trunc, tumor_data_array[xx,yy,], xlab="min", ylim=c(0, max(tumor_data_array[xx,yy,])), ylab="micro-mol/L", main="CC (nonparametric smoother)", cex=1.5)
lines(tumor_smooth_array[xx,yy,], lwd=2)

dev.set(7)
plot(map.h[xx,yy,], xlab="sec", ylim=c(0, max(map.h[xx,yy,])), ylab=expression(sec^-1), main=paste("Impulse Response Function: x=", xx, ", y=", yy, sep=""), col="white")
lines(map.h[xx,yy,], lwd=2, col="goldenrod")
text(max(TIME_trunc)/2, max(map.h[xx,yy,])*0.93, paste(paste("AUC =", format(map.AUC[xx,yy], digits=3))))
text(max(TIME_trunc)/2, max(map.h[xx,yy,])*0.86, paste("AUC/MRT =", format(map.AUCMRT[xx,yy], digits=3), "1/min"))

dev.set(4)

if(xx < (x_min + 0.25*xdim) && yy < (y_min+0.11*ydim)){
dev.set(6)
plot(TIME_trunc, TUMOR.median_corr_shifted, xlab="min", ylim=c(0, max(TUMOR.median_corr_shifted)), ylab="micro-mol/L", main="MEDIAN CC PROFILE", cex=1.5)
lines(tumor_smooth_median, lwd=2)

dev.set(7)
plot(h.median, xlab="sec", ylim=c(0, max(h.median)), ylab=expression(sec^-1), main="MEDIAN IRF", col="white")
lines(h.median, lwd=2, col="goldenrod")
text(max(TIME_trunc)/2, max(h.median)*0.93, paste(paste("AUC =", format(AUC.median, digits=3))))
text(max(TIME_trunc)/2, max(h.median)*0.86, paste("AUC/MRT =", format(AUCMRT.median, digits=3), "1/min"))

dev.set(4)
}
}
}






if(mat.original==TRUE){
cat("writing results to file...", "\n")
ptm <- proc.time()[3]


##############################################################
####SAVE RESULTS IN A MATLAB FILE#############################
##############################################################
proc.time.total <- format((proc.time()[3]-ptm_total)/60, digits=2)

###CONVERT AUCMRT UNITS FROM 1/SEC TO 1/MIN##############
#####move this to before truncation correction
#map.AUCMRT <- map.AUCMRT*60

names(proc.time.total) <- "minutes"
filename3 <- sub(".RData", "", filename3)

args <- list(file, slice, vp, border, maxCt, parameter.plot, cutoff.map, range.map, export.matlab, batch.mode, alpha.AIF, correct.trunc, kep.nom)

names(args) <- c("file", "slice", "vp", "border", "maxCt", "parameterplot", "cutoffmap", "rangemap", "exportmatlab", "batchmode", "alphaAIF", "correcttrunc", "kepnom")

plotParams <- list(nx, ny, nt, x_min, x_max, y_min, y_max, MAP_ul)
names(plotParams) <- c("nx", "ny", "nt", "xmin", "xmax", "ymin", "ymax", "MAPul")

dcemri.data <- list(vector.times, mask.roi, map_cc_slice, AIF, map.AUC, map.AUCMRT, map.h, TIME_trunc, tumor_data_array, tumor_smooth_array, ARTERY_trunc, u_smooth, TUMOR.median_corr_shifted, tumor_smooth_median, h.median, AUC.median, AUCMRT.median, args, plotParams, proc.time.total, filename3, DAT.version)

names(dcemri.data) <- c("vectorTimes", "maskROI", "mapCC", "vectorAIF", "mapAUC", "mapAUCMRT", "mapIRF", "vectorTimesTrunc", "mapCCtransformed", "mapCCsmoothed", "vectorAIFtrunc", "vectorAIFsmoothed", "CCmedianTransformed", "CCmedianSmoothed", "IRFmedian", "AUCmedian", "AUCMRTmedian", "args", "plotParams", "procTime", "filenameS", "DATversion")

save(dcemri.data, file=paste(strsplit(filename3,split=".mat")[[1]], ".RData", sep=""))

if(export.matlab==TRUE)
  writeMat(filename3, mapCC=map_cc_slice, vectorTimes=vector.times, vectorAIF=AIF, maskROI=mask.roi, mapAUC=map.AUC, mapAUCMRT=map.AUCMRT, mapCCtransformed=tumor_data_array, mapCCsmoothed=tumor_smooth_array, vectorAIFsmoothed=u_smooth, vectorAIFtrunc=ARTERY_trunc, vectorTimesTrunc=TIME_trunc, mapIRF=map.h, args=args, procTime=proc.time.total, plotParams=plotParams, DATversion=DAT.version, CCmedianSmoothed=tumor_smooth_median, IRFmedian=h.median, AUCmedian=AUC.median, AUCMRTmedian=AUCMRT.median, CCmedianTransformed=TUMOR.median_corr_shifted)

cat("..deconvolution results saved as", paste(strsplit(filename3,split=".mat")[[1]], ".RData", sep=""), "\n")

if(export.matlab==TRUE)
  cat("..deconvolution results saved as", filename3, "\n")

return(filename3)

}
}

if(file=="nodata"){
  cat("###################################################################", "\n")
  cat("###################################################################", "\n")
  cat("no data file specified...", "\n")
  cat("...loading simulated data set...", "\n")

  ## Lazy Load instead
  ## data(DAT.simData)

  cat("...writing simulated data to file...", "\n")
  dcemri.data <- DAT.simData
  save(dcemri.data, file="simulatedDCEMRI.RData")
  rm(dcemri.data)
  cat("..simulated data saved as simulatedDCEMRI.RData.", "\n", sep="")
  cat("###################################################################", "\n")
  cat("###################################################################", "\n")
  file <- "simulatedDCEMRI.RData"
  cutoff.map <- 0.95
  range.map <- 1.05
  slice <- 1
}

load(file)
data <- dcemri.data
rm(dcemri.data)
if(length(names(data)) < 8)
  mat.original <- TRUE
if(length(names(data)) >= 8)
  mat.original <- FALSE

filename3 <- DATrun(file, slice, vp, border, maxCt, parameter.plot, cutoff.map, range.map, export.matlab, batch.mode, alpha.AIF, correct.trunc)
if(batch.mode==FALSE & mat.original==TRUE){
  filetorun <- paste(strsplit(filename3,split=".mat")[[1]], ".RData", sep="")
  DATrun(filetorun, slice, vp, border, maxCt, parameter.plot, cutoff.map, range.map, export.matlab, batch.mode, alpha.AIF, correct.trunc)
}
}

