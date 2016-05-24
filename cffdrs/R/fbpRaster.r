fbpRaster <- function(input, output = "Primary", select=NULL, m=NULL, cores=1){
  #############################################################################
  # Description:
  #   An internal function used to setup the calculation of the Fire Behavior 
  #   Prediction (FBP) system over Rasters.This function moves the logic of the 
  #   FBP system equations into FBPcalc.R and sets up the use of that function 
  #   here.
  #
  #
  # Args:
  #   input:  Data frame of required and optional information needed to 
  #           calculate FBP function. View the arguments section of the fbp 
  #           manual (fbp.Rd) under "input" for the full listing of the 
  #           required and optional inputs.
  #   output: What fbp outputs to return to the user. Options are "Primary", 
  #           "Secondary" and "All".
  #   select: Select the outputs by name Ex: (select=c("HFI", "TFC", "ROS"))
  #   m:      Optimal number of pixels at each iteration of computation.
  #   cores:  Number of cores to use to parallelize this function.
  #
  # Returns:  
  #   output: Either Primary, Secondary, or all FBP outputs in a raster stack
  #
  #############################################################################
  
  #Quite often users will have a data frame called "input" already attached
  #  to the workspace. To mitigate this, we remove that if it exists, and warn
  #  the user of this case. This is also dont in FBPcalc, but we require use
  #  of this variable here before it gets to FBPCalc
  if (!is.na(charmatch("input", search()))) {
    warning("Attached dataset 'input' is being detached to use fbp() function.")
    detach(input)
  }
  #split up large rasters to allow calculation. This will be used in the
  #  parallel methods
  if (is.null(m)){
    m <- ifelse(ncell(input) > 500000, 3000, 1000)
  }
  #Setup correct output names
  allNames <- c("CFB","CFC","FD","HFI","RAZ","ROS","SFC","TFC","BE","SF","ISI",
                "FFMC", "FMC","D0", "RSO","CSI","FROS","BROS","HROSt","FROSt",
                "BROSt","FCFB", "BCFB","FFI","BFI", "FTFC","BTFC","TI","FTI",
                "BTI","LB","LBt","WSV", "DH","DB","DF","TROS","TROSt", "TCFB",
                "TFI","TTFC","TTI")
  primaryNames <- allNames[1:8]
  secondaryNames <- allNames[9:length(allNames)]
  #If outputs are specified, then check if they exist and stop with an error
  #  if not.
  if (!is.null(select)){
    select <- toupper(select)
    select <- select[!duplicated(select)]
    if(output == "SECONDARY" | output == "S"){
      if (!sort(select %in% secondaryNames)[1]){
        stop("Selected variables are not in the outputs")}
    }
    if (output == "PRIMARY" | output == "P"){
      if (!sort(select %in% primaryNames)[1]){
        stop("Selected variables are not in the outputs")} 
    }
    if (output == "ALL" | output == "A"){
      if (!sort(select %in% allNames)[1]){
        stop("Selected variables are not in the outputs")} 
    }
  }
  names(input) <- toupper(names(input))
  output <- toupper(output)
  if("LAT" %in% names(input)){
    #register a sequential parallel backend
    registerDoSEQ()
    #Get the specified raster cell values
    r <- getValuesBlock_stackfix(input, nrows=nrow(input))
    #convert to data.frame
    r <- as.data.frame(r)
    names(r) <- names(input)
  }else{
    #Convert the raster to points and insert into a data.frame
    r <- as.data.frame(rasterToPoints(input))
    #Rename the latitude field
    names(r)[names(r) == "y"] <- "LAT"
    #Check for valid latitude
    if (max(r$LAT) > 90 | min(r$LAT) < -90){
      warning("Input projection is not in lat/long, consider re-projection or 
              include LAT as input")
    }
  }
  #Unique IDs
  r$ID <- 1:nrow(r)
  #merge fuel codes with integer values
  fuelCross <- data.frame(FUELTYPE0 = sort(c(paste("C", 1:7, sep="-"),
                                             "D-1",
                                             paste("M", 1:4, sep="-"),
                                             paste("S", 1:3, sep="-"),
                                             "O-1a", "O-1b", "WA", "NF")),
                          code=1:19)
  r <- merge(r, fuelCross, by.x="FUELTYPE", by.y="code", all.x=TRUE, all.y=FALSE)

  r$FUELTYPE <- NULL
  names(r)[names(r) == "FUELTYPE0"] <- "FUELTYPE"
  r <- r[with(r, order(ID)), ]
  names(r)[names(r) == "x"] <- "LONG"
  #Calculate FBP through the fbp() function
  FBP <- fbp(r, output = output, m = m, cores = cores)
  #If secondary output selected then we need to reassign character
  #  represenation of Fire Type S/I/C to a numeric value 1/2/3
  if (!(output == "SECONDARY" | output == "S")){
    FBP$FD <- ifelse(FBP$FD == "I", 2, FBP$FD)
    FBP$FD <- ifelse(FBP$FD == "C", 3, FBP$FD)
    FBP$FD <- ifelse(FBP$FD == "S", 1, FBP$FD)
    FBP$FD <- as.numeric(FBP$FD)
  }
  #If caller specifies select outputs, then create a raster stack that contains
  #  only those outputs
  if (!is.null(select)){
    out <- out0 <- input[[1]]
    values(out) <- FBP[, select[1]]
    if (length(select) > 1){
      for (i in 2:length(select)){
        values(out0) <- FBP[,select[i]]
        out <- stack(out, out0)
      }
    }
    names(out)<-select 
  #If caller specified Primary outputs, then create raster stack that contains
  #  only primary outputs
  }else if (output == "PRIMARY" | output == "P") {
    message("FD = 1,2,3 representing Surface (S),Intermittent (I), and Crown (C) fire")
    out <- out0 <- input[[1]]
    values(out) <- FBP[,primaryNames[1]]
    for (i in 2:length(primaryNames)){
      values(out0) <- FBP[, primaryNames[i]]
      out <- stack(out,out0)
    }
    names(out)<-primaryNames
  #If caller specified Secondary outputs, then create raster stack that contains
  #  only secondary outputs
  }else if(output == "SECONDARY" | output == "S") {
    out <- out0 <- input[[1]]
    values(out) <- FBP[, secondaryNames[1]]
    for (i in 2:length(secondaryNames)){
      values(out0) <- FBP[, secondaryNames[i]]
      out <- stack(out, out0)
    }
    names(out)<-secondaryNames
  #If caller specified All outputs, then create a raster stack that contains
  #  both primary and secondary outputs
  }else if(output == "ALL" | output == "A") {
    message("FD = 1,2,3 representing Surface (S),Intermittent (I), and Crown (C) fire")
    out <- out0 <- input[[1]]
    values(out) <- FBP[, allNames[1]]
    for (i in 2:length(allNames)){
      values(out0) <- FBP[, allNames[i]]
      out <- stack(out, out0)
    }
    names(out) <- allNames
  }
  #return the raster stack to the caller
  return(out)
}
