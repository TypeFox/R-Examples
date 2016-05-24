QualityCheck <- 
function(Data, Product, Band, NoDataFill, QualityBand, QualityScores, QualityThreshold)
{
  ##### Define what valid ranges for quality bands should be for different products:
  # 1=lower range     2=upper range     3=no fill value
  QA_RANGE <- data.frame(
    MOD09A1 = c(0, 4294966531, ""),            # Surface reflectance bands 0-1 bits
    MYD09A1 = c(0, 4294966531, ""),            # Surface reflectance bands 0-1 bits
    MOD11A2 = c(0, 255, ""),                   # Land surface temperature and emissivity 0-1 bits
    MYD11A2 = c(0, 255, ""),                   # Land surface temperature and emissivity 0-1 bits
    MCD12Q1 = c(0, 254, 255),                  # Land cover types 0-1 bits
    MOD13Q1 = c(0, 3, -1),                     # Vegetation indices 0-1 bits
    MYD13Q1 = c(0, 3, -1),                     # Vegetation indices 0-1 bits
    MOD15A2 = c(0, 254, 255),                  # LAI - FPAR 0 bit
    MYD15A2 = c(0, 254, 255),                  # LAI - FPAR 0 bit
    MOD17A2 = c(0, 254, 255),                  # GPP 0 bit
    MOD17A3 = c(0, 100, ""),                   # GPP funny one - evaluate separately.
    MCD43A4 = c(0, 4294967294, 4294967295)     # BRDF albedo band quality, taken from MCD43A2, for reflectance data.
  )                  
  # Land cover dynamics products are available for download but not for quality checking with this function.
  
  # Check the product input corresponds to one with useable quality information
  if(!any(names(QA_RANGE) == Product)) stop(paste("QualityCheck cannot be used for (", Product, ") product.", sep = ""))
  
  # If dataframes, coerce to matrices.
  if(is.data.frame(Data)) Data <- as.matrix(Data)
  if(is.data.frame(QualityScores)) QualityScores <- as.matrix(QualityScores)  
  
  # Check that Data and QualityScores have matching length and, if a matrix, dimensions .
  if(is.matrix(Data) | is.matrix(QualityScores)){
    if(!(is.matrix(Data) & is.matrix(QualityScores))) stop("Data and QualityScores do not have matching dimensions.")
    if(!all(dim(Data) == dim(QualityScores))) stop("Data and QualityScores do not have matching dimensions.")
  } else {
    if(length(Data) != length(QualityScores)) stop("Data and QualityScores must have matching lengths.")
  }
  
  # Check the QualityScores input are within the correct range for the product requested.
  if(any((QualityScores < QA_RANGE[1,Product] | QualityScores > QA_RANGE[2,Product]) & 
          QualityScores != QA_RANGE[3,Product])){
    stop(paste("Some QualityScores are outside of the range of valid quality control data for the product requested. For 
         this product, the valid QualityScore range is", QA_RANGE[1,Product], "and", QA_RANGE[2,Product], ".", sep = " "))
  }
  
  # Quality Threshold should be one integer.
  if(length(QualityThreshold) != 1) stop("QualityThreshold input must be one integer.")
  if(!is.numeric(QualityThreshold)) stop("QualityThreshold should be numeric class. One integer.")
  if(abs(QualityThreshold - round(QualityThreshold)) > .Machine$double.eps^0.5){
    stop("QualityThreshold input must be one integer.")
  }
  
  # NoDataFill should be one integer.
  if(length(NoDataFill) != 1) stop("NoDataFill input must be one integer.")
  if(!is.numeric(NoDataFill)) stop("NoDataFill should be numeric class. One integer.")
  if(abs(NoDataFill - round(NoDataFill)) > .Machine$double.eps^0.5) stop("NoDataFill input must be one integer.")
  
  # Check QualityThreshold is within 0-3 for all except LAI and GPP bands, which must be within 0-1, and MOD17A3 (0-100).
  lai.gpp.prods <- c("MOD15A2", "MYD15A2", "MOD17A2")
  if(any(lai.gpp.prods == Product)){
    if(QualityThreshold != 0 & QualityThreshold != 1){
      stop("QualityThreshold should be either 0 or 1 for this product; 0 is good quality, 1 is other.")
    }
  } else if(Product != "MOD17A3"){
    if(QualityThreshold < 0 | QualityThreshold > 3){
      stop("QualityThreshold should be between 0-3 for this product; 
           0 = high quality, 1 = good but marginal quality, 2 = cloudy/poor quality, 3 = poor quality for other reasons.")
    }
  }
  #####
  
  # MOD17A3 is an exception, so deal with this first, and then the rest,
  if(Product == "MOD17A3"){
    if(QualityThreshold < 0 | QualityThreshold > 100) stop("QualityThreshold should be between 0-100 for this product.")
    
    NOFILLRANGE <- c(250, 255)
    Data <- ifelse(Data > NOFILLRANGE[1] & Data < NOFILLRANGE[2], Data, NA)
    Data <- ifelse(QualityScores <= QualityThreshold, Data, NA)
  } else {
    # Convert decimal QualityScores values into binary.
    decimal.set <- QualityScores
    
    if(max(QualityScores) == 0) num.binary.digits <- 1
    if(max(QualityScores) != 0) num.binary.digits <- floor(log(max(QualityScores), base = 2)) + 1
    
    binary.set<- matrix(nrow = length(QualityScores), ncol = num.binary.digits)
    
    for(n in 1:num.binary.digits){
      binary.set[ ,(num.binary.digits - n) + 1] <- decimal.set %% 2
      decimal.set <- decimal.set %/% 2
    }
    
    quality.binary <- apply(binary.set, 1, function(x) paste(x, collapse = ""))
    
    # Deal with the quality data in binary, according to the product input.
    if(any(lai.gpp.prods == Product)){
      qa.binary <- as.numeric(substr(quality.binary, nchar(quality.binary), nchar(quality.binary)))
      Data <- ifelse(Data != NoDataFill & qa.binary <= QualityThreshold, Data, NA)
    } else {
      # Create an ifelse here so if MCD43A4, snip the relevant binary string for Band. Otherwise, carry on.
      if(Product == "MCD43A4"){
        band.num <- as.numeric(substr(Band, nchar(Band), nchar(Band)))
        if(any(is.na(band.num))) stop("Band input is not one of the reflectance bands (1-7) from MCD43A4.")
        if(any(1 < band.num & band.num > 7)) stop("Band input is not one of the reflectance bands (1-7) from MCD43A4.")
        
        # Select the section of binary code relevant to Band.
        qa.binary <- substr(quality.binary, (nchar(quality.binary) - (((band.num - 1) * 2) + 2)), 
                                            (nchar(quality.binary) - ((band.num - 1) * 2)))

        qa.int <- numeric(length(qa.binary))
        qa.int[qa.binary == "000"] <- 0
        qa.int[qa.binary == "001"] <- 1
        qa.int[qa.binary == "010"] <- 2
        qa.int[qa.binary == "011"] <- 3
        qa.int[qa.binary == "100"] <- 4
        
        Data <- ifelse(Data != NoDataFill & qa.int <= QualityThreshold, Data, NA)
      } else {
        qa.binary <- substr(quality.binary, nchar(quality.binary) - 1, nchar(quality.binary))

        qa.int <- numeric(length(qa.binary))   
        qa.int[qa.binary == "00"] <- 0
        qa.int[qa.binary == "01"] <- 1
        qa.int[qa.binary == "10"] <- 2
        qa.int[qa.binary == "11"] <- 3
        
        Data <- ifelse(Data != NoDataFill & qa.int <= QualityThreshold, Data, NA)
      }
    }
  }
  return(Data)
}