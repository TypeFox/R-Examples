#' Run a batch of molaR analyses
#'
#' A function which automats molaR analyses. User simply sets up the functions they
#' want run and can leave the computer to do the rest. 
#'
#' @param pathname The path to the file containing all the PLY surfaces to be
#' analyzed. Defaults to the working directory
#' @param filename Name for the output csv file. Users must include .csv suffix
#' @param DNE logical indicating whether or not to perform DNE calculation Defaults
#' to true
#' @param RFI logical indicating whether or not to perform RFI calculation Defaults to
#' true
#' @param OPCr logical indicating whether or not to perform OPCr calculation
#' Defaults to true
#' @param OPC logical indicating whether or not to perform OPC calculation Defaults
#' to false
#' @param Details logical indicating whether or not to save the details of the RFI and
#' OPCr calculations
#' @param DNE_outliers the percentile at which outliers will be excluded is passed to
#' the DNE function, defaults to 0.1
#' @param DNE_BoundaryDiscard is a logical indicating how
#' to handle the exclusion of the faces on the edge of the surface, defaults to excluding
#' faces which have a leg on the boundary. 
#' @param RFI_alpha the size of the alpha passed to RFI
#' function, defaults to 0.01
#' @param OPCr_steps the number of steps the OPCr function should take, is passed to
#' the OPCr function. Defaults to 8
#' @param OPCr_stepSize the size of each rotation. Passed to the OPCr function.
#' Defaults to 5.626 degrees
#' @param OPCr_minimum_faces sets the lower boundary for number of faces a patch
#' must have for inclusion in total count. Defaults to 3 or more.
#' @param OPCr_minimum_area sets the lower boundary for percentage of the surface
#' area a patch must make up for inclusion in the total patch count. Cannot be used with
#' minimum_faces on. Defaults to zero
#' @param OPC_rotation amount of rotation to apply during OPC calculation. Defaults
#' to zero
#' @param OPC_minimum_faces minimum number of faces a patch must contain to be
#' counted in the OPC function. Defaults to 3.
#' @param OPC_minimum_area minimum percentage of the surface area a patch must
#' make up to be counted in the OPC function. Defaults to off
#' @param Parameters defaults to off. When engaged a list of all the parameters used
#' during molaR analysis will be appended to the output file. 
#'
#' @details This function allows a user to set the analyses from molaR they want to run,
#' along with the specific parameters for each function and have a whole batch of PLY
#' files analyzed and saved to a csv file. Function will perform analyses on all PLY files
#' in the working directory or user can specify a file path.
#'
#' @importFrom
#' Rvcg vcgPlyRead
#'
#' @export
#' molaR_Batch


molaR_Batch <- function (pathname = getwd(), filename = "molaR_Batch.csv", DNE = TRUE, 
          RFI = TRUE, OPCr = TRUE, OPC = FALSE, Details = FALSE, DNE_outliers = 0.1, 
          DNE_BoundaryDiscard = "Leg", RFI_alpha = 0.01, OPCr_steps = 8, 
          OPCr_stepSize = 5.626, OPCr_minimum_faces = 3, OPCr_minimum_area = 0, 
          OPC_rotation = 0, OPC_minimum_faces = 3, OPC_minimum_area = 0, Parameters = FALSE) 
{
  if (DNE == FALSE && RFI == FALSE && OPCr == FALSE && OPC == FALSE) {
    stop("No metrics were selected to run")
  }
  fileNames <- dir(pathname, pattern = "*.ply")
  if (length(fileNames) == 0) {
    stop("No PLY files in this directory")
  }
  fileNum <- length(fileNames)
  DNE_Output <- vector("numeric")
  RFI_Output <- vector("numeric")
  OPCr_Output <- vector("numeric")
  OPC_Output <- vector("numeric")
  if (Details == TRUE) {
    ThreeD_Area <- vector("numeric")
    TwoD_Area <- vector("numeric")
    RotDegrees <- seq(0, (OPCr_steps * OPCr_stepSize) - OPCr_stepSize, 
                      OPCr_stepSize)
    Rotations <- matrix(nrow = length(fileNames), ncol = length(RotDegrees))
    colnames(Rotations) <- paste(as.character(RotDegrees), 
                                 "deg.")
  }
  for (i in 1:fileNum) {
    invisible(capture.output(Specimen <- molaR_Clean(vcgPlyRead(file.path(pathname, 
                                                                          fileNames[i])))))
    if (DNE == TRUE) {
      invisible(capture.output({
        DNE_Specimen <- try(DNE(Specimen, outliers = DNE_outliers, 
                                BoundaryDiscard = DNE_BoundaryDiscard))
      }))
      if (is.character(DNE_Specimen)) {
        DNE_Output = c(DNE_Output, DNE_Specimen)
      }
      else {
        DNE_Result <- DNE_Specimen$Surface_DNE
        DNE_Output <- c(DNE_Output, DNE_Result)
      }
    }
    if (RFI == TRUE) {
      invisible(capture.output({
        RFI_Specimen <- try(RFI(Specimen, alpha = RFI_alpha))
      }))
      if (is.character(RFI_Specimen)) {
        RFI_Output <- c(RFI_Output, RFI_Specimen)
        if (Details == TRUE) {
          ThreeD_Result <- NA
          TwoD_Result <- NA
          ThreeD_Area <- c(ThreeD_Area, ThreeD_Result)
          TwoD_Area <- c(TwoD_Area, TwoD_Result)
        }
      }
      else {
        RFI_Result <- RFI_Specimen$Surface_RFI
        RFI_Output <- c(RFI_Output, RFI_Result)
        if (Details == TRUE) {
          ThreeD_Result <- RFI_Specimen$Three_D_Area
          TwoD_Result <- RFI_Specimen$Two_D_Area
          ThreeD_Area <- c(ThreeD_Area, ThreeD_Result)
          TwoD_Area <- c(TwoD_Area, TwoD_Result)
        }
      }
    }
    if (OPCr == TRUE) {
      if (OPCr_steps != round(OPCr_steps)) {
        stop("Please enter integer value for number of OPCr steps")
      }
      invisible(capture.output({
        OPCr_Specimen <- try(OPCr(Specimen, Steps = OPCr_steps, 
                                  stepSize = OPCr_stepSize, minimum_faces = OPCr_minimum_faces, 
                                  minimum_area = OPCr_minimum_area))
      }))
      if (is.character(OPCr_Specimen)) {
        OPCr_Output = c(OPCr_Output, OPCr_Specimen)
        if (Details == TRUE) {
          for (j in 1:OPCr_steps) {
            Rotations[i, j] <- NA
          }
        }
      }
      else {
        OPCr_Result <- OPCr_Specimen$OPCR
        OPCr_Output <- c(OPCr_Output, OPCr_Result)
        if (Details == TRUE) {
          for (j in 1:OPCr_steps) {
            Rotations[i, j] <- OPCr_Specimen$Each_Run[[j, 
                                                       2]]
          }
        }
      }
    }
    if (OPC == TRUE) {
      invisible(capture.output({
        OPC_Specimen <- OPC(Specimen, rotation = OPC_rotation, 
                            minimum_faces = OPC_minimum_faces, minimum_area = OPC_minimum_area)
      }))
      OPC_Result <- OPC_Specimen$Patch_Count$`total patches`
      OPC_Output <- c(OPC_Output, OPC_Result)
    }
    cat("\nProcessed", i, "of", fileNum, "total PLY files in directory")
  }
  Output <- data.frame(Files = fileNames)
  if (DNE == TRUE) {
    Output <- cbind(Output, DNE = DNE_Output)
  }
  if (RFI == TRUE) {
    Output <- cbind(Output, RFI = RFI_Output)
    if (Details == TRUE) {
      Output <- cbind(Output, `3D_Area` = ThreeD_Area, 
                      `2D_Area` = TwoD_Area)
    }
  }
  if (OPCr == TRUE) {
    Output <- cbind(Output, OPCR = OPCr_Output)
    if (Details == TRUE) {
      Output <- cbind(Output, Rotations)
    }
  }
  if (OPC == TRUE) {
    Output <- cbind(Output, OPC = OPC_Output)
  }
  print(Output)
  if (Parameters == TRUE){
    paramList <- matrix(data="", nrow=2, ncol=ncol(Output))
    colnames(paramList) <- colnames(Output)
    paramList[2,1] <- "Parameter"
    paramList[2,2] <- "Value"
    Blanks <- ncol(paramList)-2
    if (DNE == TRUE) {
      newParam1 <- c("DNE_outliers", DNE_outliers, rep("", Blanks))
      newParam2 <- c("DNE_BoundaryDiscard", DNE_BoundaryDiscard, rep("", Blanks))
      paramList <- rbind(paramList, newParam1, newParam2)
    }
    if (RFI == TRUE) {
      newParam3 <- c("RFI_alpha", RFI_alpha, rep("", Blanks))
      paramList <- rbind(paramList, newParam3)
    }
    if (OPCr == TRUE) {
      newParam4 <- c("OPCr_steps", OPCr_steps, rep("", Blanks))
      newParam5 <- c("OPCr_stepSize", OPCr_stepSize, rep("", Blanks))
      newParam6 <- c("OPCr_minimum_faces", OPCr_minimum_faces, rep("", Blanks))
      newParam7 <- c("OPCr_minimum_area", OPCr_minimum_area, rep("", Blanks))
      paramList <- rbind(paramList, newParam4, newParam5, newParam6, newParam7)
    }
    if (OPC ==TRUE) {
      newParam8 <- c("OPC_rotation", OPC_rotation, rep("", Blanks))
      newParam9 <- c("OPC_minimum_faces", OPC_minimum_faces, rep("", Blanks))
      newParam10 <- c("OPC_minimum_area", OPC_minimum_area, rep("", Blanks))
      paramList <- rbind(paramList, newParam8, newParam9, newParam10)
    }
    Output <- rbind(Output, paramList)
  }
  write.csv(Output, file = file.path(pathname, filename), row.names = FALSE)
  cat("Results saved to directory.")
}