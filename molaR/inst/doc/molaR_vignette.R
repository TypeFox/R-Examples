## ----setup, message = FALSE----------------------------------------------
library(molaR)
summary(ex_tooth1)

## ----webgl_setup, echo=FALSE---------------------------------------------
library(knitr)
library(rgl)

## ----DNE_basic, webgl = TRUE---------------------------------------------
DNE1 = DNE(ex_tooth1)
DNE3d(DNE1) #you will need to zoom out to see the whole tooth in an html graphics window

## ----DNE_color_scale, webgl = TRUE---------------------------------------
DNE2 = DNE(ex_tooth2)
DNE3d(DNE2, setRange = c(0, 1.3))

## ----DNE_color_scale2, webgl = TRUE--------------------------------------
DNE3d(DNE1, setRange = c(0, 1.3))

## ----DNE_object----------------------------------------------------------
head(DNE1$Edge_Values)
head(DNE1$Outliers)

## ----RFI_basic-----------------------------------------------------------
RFI1 = RFI(ex_tooth1, alpha=0.5)

## ----RFI_plot, webgl = TRUE----------------------------------------------
RFI3d(RFI1)

## ----OPC_basic, webgl = TRUE---------------------------------------------
OPC1 = OPC(ex_tooth1)
OPC3d(OPC1)

## ----OPC_patch_count-----------------------------------------------------
OPC2 = OPC(ex_tooth1, minimum_faces = 20)
OPC3 = OPC(ex_tooth1, minimum_area = 0.01)

## ----OPCr----------------------------------------------------------------
OPCr1 = OPCr(ex_tooth1)
OPCr2 = OPCr(ex_tooth1, Steps = 5, stepSize = 9, minimum_faces = 2) #the minimum_faces and minimum_area parameters are passed to each iteration of OPC

## ----OPCr_structure------------------------------------------------------
OPCr1$Each_Run
OPCr2$Each_Run

