## ----, echo = FALSE, message = FALSE-------------------------------------
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)

## ----, message = FALSE---------------------------------------------------
library(hazus)
library(reshape2)
library(ggplot2)

## ------------------------------------------------------------------------
data(haz_fl_dept) # depth-based DFs
data(haz_fl_velo) # velocity and depth-based DFs
data(haz_fl_agri) # agriculture DFs
data(haz_fl_bridge) # DFs for bridges
data(haz_fl_depr) # depreciation functions
data(haz_fl_occ) # occupancy description table

## ------------------------------------------------------------------------
fl_dept <- extract_hazus_functions(long_format = FALSE)
dim(fl_dept)

## ------------------------------------------------------------------------
fl_dept <- extract_hazus_functions()
dim(fl_dept)
head(fl_dept)

## ------------------------------------------------------------------------
head(haz_fl_occ)
levels(as.factor(haz_fl_occ$Occupancy))
table(haz_fl_occ$Occupy_Class)

## ------------------------------------------------------------------------
fl_velo <- extract_hazus_functions(func_type = "velocity")
str(fl_velo)

## ------------------------------------------------------------------------
fl_agri <- extract_hazus_functions(func_type = "ag")
str(fl_agri)

## ------------------------------------------------------------------------
fl_bridge <- extract_hazus_functions(func_type = "bridge")
str(fl_bridge)

## ------------------------------------------------------------------------
fl_depr <- extract_hazus_functions(func_type = "deprec")
str(fl_depr)

## ------------------------------------------------------------------------
gfx_data <- subset(fl_dept, grepl("one floor", Description) & Cover_Class == "Bldg")

# clean up description
gfx_data$Description <- paste(gfx_data$DmgFnId, gfx_data$Description)

gfx_line <- ggplot(data = gfx_data, aes(x = depth, y = damage))
gfx_line <- gfx_line + geom_line(aes(colour = Description))
gfx_line <- gfx_line + ylab("Damage (%)") + xlab("Flood Depth (ft)")

print(gfx_line)

## ------------------------------------------------------------------------
gfx_line <- ggplot(data = fl_velo, aes(x = vel, y = dep))
gfx_line <- gfx_line + geom_line(aes(colour = num_story))
gfx_line <- gfx_line + facet_wrap(~struc_type, scales = "fixed")
gfx_line <- gfx_line + ylab("Flood Depth (ft)") + xlab("Flood Velocity (ft/s)")
print(gfx_line)

## ------------------------------------------------------------------------
gfx_data <- subset(fl_agri, loss_type == "PctCropLoss" & Crop %in% 
                     c("Tomato", "Cotton", "Wheat"))

gfx_line <- ggplot(data = gfx_data, aes(x = JulianDay, y = damage))
gfx_line <- gfx_line + geom_line(aes(colour = Crop))
gfx_line <- gfx_line + ylab("Damage (fraction)") + xlab("Day of Year")
print(gfx_line)

## ------------------------------------------------------------------------
gfx_data <- fl_bridge
gfx_data$Description <- paste(gfx_data$Occupancy, gfx_data$Description)

gfx_line <- ggplot(data = gfx_data, aes(x = ret_period, y = damage))
gfx_line <- gfx_line + geom_line(aes(colour = Description))
gfx_line <- gfx_line + ylab("Damage (%)") + xlab("Flood Return Period (years)")
print(gfx_line)

## ------------------------------------------------------------------------
gfx_data <- fl_depr[grepl("1", fl_depr$Occupancy), ]

gfx_line <- ggplot(data = gfx_data, aes(x = Age, y = deprec))
gfx_line <- gfx_line + geom_line(aes(colour = Occupancy))
gfx_line <- gfx_line + ylab("Depreciation (%)") + xlab("Age (years)")
print(gfx_line)

