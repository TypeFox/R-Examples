##  SOIL CLASSIFICATION
##  Kyle Elmy and Jim Kaklamanos
##  1 December 2015




########################################################################################
##  1. AASHTO SOIL CLASSIFICATION
##
##  Description:  This set of functions classifies soil using the American Association
##                of State Highway and Transportation Officials (AASHTO)
##                Soil Classification System

##------------------------------------------------------------------------------------------


##------------------------------------------------------------------------------------------
##  1a. AASHTO GROUP CLASSIFICATION

##  Description:   Inorganic soil classification using the AASHTO method
##
##  Output:   AASHTO group classification
##
##  Inputs:   Grain-size parameters:
##            o  sieve = vector of sieve numbers (according to ASTM D422)
##                       in grain-size distribution
##            o  size = vector of grain sizes (in or mm) in distribution
##            o  percent = vector of percent passing in grain-size distribution
##            o  metric = logical variable for grain-size distribution:
##                        TRUE for metric units (mm), FALSE for English units (in);
##                        only required if "size" is supplied
##            o  p10, p40, p200 = percent passing #10 sieve, #40 sieve, and #200 sieve
##
##            Plasticity parameters:
##            o  LL = liquid limit (percent)
##            o  PL = plastic limit (percent)
##            o  PI = plasticity index (percent)
##            o  NP = logical variable indicating whether the soil is nonplastic (TRUE or FALSE)
##
##  Notes: o  For the grain-size data, one of three following pieces of data must be specified:
##            1.  Sieve numbers (sieve); and percent passing
##            2.  Grain sizes (size); and percent passing
##            3.  p10, p40, p200
##         o  If sieve data are specified, for sieves larger than the no. 4 sieve,
##            the user should specify the sieve size in inches 
##            (e.g., 3/8, 3/4, 1, 1.5, 2, 3, etc.)
##         o  If sieve or size data are specified, this function assumes that the
##            no. 10, 40, and 200 sieves have been used.
##         o  For the plasticity data, one of three following pieces of data must be specified:
##            1.  LL and PL
##            2.  LL and PI
##            3.  NP (if the soil is nonplastic)
##         o  To obtain the group index using the AASHTO method, use the GI function.

AASHTO <- function(sieve = NA, size = NA, percent = NA, metric = NA, p10 = NA, p40 = NA, p200 = NA,
                   LL = NA, PL = NA, PI = NA, NP = NA){

  ##  ------------------------------------------------------------
  ##  Plasticity parameters
  
  ##  Calculate PI
  if(is.na(PI) == TRUE && is.na(LL) == FALSE && is.na(PL) == FALSE){
    PI <- LL - PL
  }

  ##  Non-plastic soils
  if(is.na(NP) == FALSE){
    if(NP == TRUE){
      PI <- 0
      LL <- 0
      index <- 0
    }
  }

  ##  ------------------------------------------------------------
  ##  Grain-size parameters
  
  ##  Obtain grain sizes from sieves, if sieve numbers are provided
  if(all(is.na(sieve) == FALSE)){
    size <- size.from.sieve(sieve = sieve, metric = TRUE)
    metric <- TRUE
  }

  ##  Obtain p10, p40, p200 if size or sieve numbers are provided
  if(all(is.na(size) == FALSE)){
    if(metric == TRUE){
      p10 <- percent[match(2, round(size, 1))]
      p40 <- percent[match(0.425, round(size, 3))] 
      p200 <- percent[match(0.075, round(size, 3))]
    } else{
      if(metric == FALSE){
        p10 <- percent[match(0.079, round(size, 3))]
        p40 <- percent[match(0.017, round(size, 3))] 
        p200 <- percent[match(0.003, round(size, 3))]
      }
    }
  }
  
  ##  ------------------------------------------------------------
  ##  Coarse-grained soil classification

  ##  Stone fragments; gravel and sand
  if(p200 <= 35){
    if(p10 <= 50 && p40 <= 30 && p200 <= 15 && PI <= 6){
      class <- "A-1-a"
    } else if(p40 <= 50){
      if(p200 <= 25 && PI <= 6){
        class <- "A-1-b"
          
      ##  Classify silty or clayey gravel & sand
      } else if(PI <= 10){
        if(LL <= 40){
          class <- "A-2-4"
        } else if(LL > 40){
          class <- "A-2-5"
        }
      } else if(PI > 10){
        if(LL <= 40){
          class <- "A-2-6"
        } else if(LL > 40){
          class <- "A-2-7"
        }
      }
      
      ## Classify fine sand
    } else if(p40 > 50){
      if(p200 <= 10 && NP == TRUE){
        class <- "A-3"
      } else{
        if(PI <= 10){
          if(LL <= 40){
            class <- "A-2-4"
          } else if(LL > 40){
            class <- "A-2-5"
          }
        } else if(PI > 10){
          if(LL <= 40){
            class <- "A-2-6"
          } else if(LL > 40){
            class <- "A-2-7"
          }
        }          
      }
    }
  
  ##  ------------------------------------------------------------
  ##  Fine-grained soil classification
    
  } else if(p200 > 35){
    
    ## Classify silty soils
    if(PI <= 10){
      if(LL <= 40){
        class <- "A-4"
      } else if(LL > 40){
        class <- "A-5"
      }
      
    ## Classify clayey soils
    } else if(PI > 10){
      if(LL <= 40){
        class <- "A-6"
      } else if(LL > 40){
        if(PI <= (LL-30)){
          class <- "A-7-5"
        } else if(PI > (LL-30)){
          class <- "A-7-6"
        }
      }
    }
  }
  
  ##  Return class
  return(class)
}
  
##  Example:
##
##  AASHTO(p10 = 51, p40 = 30, p200 = 15, LL = 40, PI = 10)




##------------------------------------------------------------------------------------------
##  1b. AASHTO GROUP INDEX
##
##  Description:  Calculate a soil's group index using the AASHTO method
##
##  Output:  AASHTO group index (rounded to the nearest integer)
##
##  Inputs:  p200 = fines content, i.e. percent passing the #200 sieve (percent)
##           LL = liquid limit (percent)
##           PL = plastic limit (percent)
##           PI = plasticity index (percent)
##
##  Notes:   Either PL or PI must be specified

GI <- function(p200, LL, PL = NA, PI = NA){
  
  ##  Calculate PI if not provided
  if(is.na(PI) == TRUE){
    PI <- LL - PL
  }

  ##  Calculate group index
  F <- p200
  GI <- ((F-35)*(0.2 + (0.005*(LL-40)))) + (0.01*(F-15)*(PI-10))
  GI <- round(GI, 0)  ##  Round to nearest integer
  GI <- max(GI, 0)    ##  Coerce negative numbers to zero

  ##  Return value of Group Index
  return(GI)
}

##  Example:
##
##  GI(F = 48, LL = 45, PI = 21)





########################################################################################
##  2. USCS SOIL CLASSIFICATION
##
##  Description:  This set of functions classifies soil using the Unified Soil
##                Classification System (USCS)

##------------------------------------------------------------------------------------------



##------------------------------------------------------------------------------------------
##  2a. USCS FINE-GRAINED SYMBOL
##
##  Output:   A soil's fine-grained USCS symbol
##
##  Inputs:   LL = Liquid Limit
##            PL = Plastic Limit
##            PI = Plasticity Index
##
##  Notes:    Either PL or PI must be provided.

USCS.fine.symbol <- function(LL, PL = NA, PI = NA){

  ##  Calculate PI
  if(is.na(PI)){
    PI <- LL - PL
  }

  ##  Make sure the point is below the 1-to-1 line
  if(PI > LL){
    stop("PI cannot be greater than LL.")
  }
  
  ##  Determine whether above or below A-line (first letter)
  if(PI >= A.line(LL)){
    Letter1 <- "C"
  } else if(PI < A.line(LL)){
    Letter1 <- "M"
  }

  ##  Determine whether high or low plasticity (second letter)
  if(LL >= 50){
    Letter2 <- "H"
  } else if(LL < 50){
    Letter2 <- "L"
  }

  ##  Obtain group symbol
  symbol <- paste(Letter1, Letter2, sep = "")
  
  ##  Deal with CL-ML case
  if(symbol == "CL"){
    if(PI <= 7){
      symbol <- "CL-ML"
    }
  }
  return(symbol)
}




##------------------------------------------------------------------------------------------
##  2b. USCS COARSE-GRAINED SYMBOL
##
##  Output:   A soil's USCS coarse-grained symbol
##
##  Inputs:   pg = percent gravel (%)
##            ps = percent sand (%)
##            pf = percent fines (%)
##            Cu = coefficient of uniformity
##            Cc = coefficient of curvature
##            LL = liquid limit
##            PL = plastic limit
##            PI = plasticity index
##
##  Notes: o  Data on the soil's gradation (Cu and Cc) are required if
##            the percent fines is less than or equal to 12%.
##         o  Data on the soil's fines [either (a) LL and PL or (b) LL and PI] are
##            required if the percent fines is greater than or equal to 5%.

USCS.coarse.symbol <- function(pg, ps, pf, Cc, Cu, PI = NA, LL, PL){

  ##  Sand vs. gravel (first letter)
  if(ps < pg){
    coarse <- "G"
  } else if(ps >= pg){
    coarse <- "S"
  }

  ##  ------------------------------
  ##  Second letter: gravel
  if(coarse == "G"){
    ##  Percent fines < 5%
    if(pf < 5){
      if(Cu >= 4){
        if(Cc >= 1 && Cc <= 3){
          symbol <- "GW"
        }
      } else{
        symbol <- "GP"
      }

      ##  5% <= Percent Fines <= 12%
    } else if(pf >= 5 && pf <= 12){
      ##  Get USCS group symbol
      fines <- USCS.fine.symbol(PI = PI, LL = LL, PL = PL)
      if(Cu >= 4){
        if(Cc >= 1 && Cc <= 3){
          if(substr(fines,1,1) == "M"){
            symbol <- "GW-GM"
          } else if(substr(fines,1,1) == "C"){
            symbol <- "GW-GC"
          }
        }
      } else{
        if(substr(fines,1,1) == "M"){
          symbol <- "GP-GM"
        } else if(substr(fines,1,1) == "C"){
          symbol <- "GP-GC"
        }
      }

      ##  Percent Fines > 12%
    } else if(pf > 12){
      fines <- USCS.fine.symbol(PI = PI, LL = LL, PL = PL)
      if(substr(fines,1,1) == "M"){
        symbol <- "GM"
      } else if(fines == "CL" || fines == "CH"){
        symbol <- "GC"
      } else if(fines == "CL-ML"){
        symbol <- "GC-GM"
      }
    }
  }
        

  ##  ------------------------------
  ##  Second letter: sand
  if(coarse == "S"){
    ##  Percent fines < 5%
    if(pf < 5){
      if(Cu >= 6){
        if(Cc >= 1 && Cc <= 3){
          symbol <- "SW"
        }
      } else{
        symbol <- "SP"
      }

      ##  5% <= Percent Fines <= 12%
    } else if(pf >= 5 && pf <= 12){
      ##  Get USCS group symbol
      fines <- USCS.fine.symbol(PI = PI, LL = LL, PL = PL)
      if(Cu >= 6){
        if(Cc >= 1 && Cc <= 3){
          if(substr(fines,1,1) == "M"){
            symbol <- "SW-SM"
          } else if(substr(fines,1,1) == "C"){
            symbol <- "SW-SC"
          }
        }
      } else{
        if(substr(fines,1,1) == "M"){
          symbol <- "SP-SM"
        } else if(substr(fines,1,1) == "C"){
          symbol <- "SP-SC"
        }
      }
      
      ##  Percent Fines > 12%
    } else if(pf > 12){
      fines <- USCS.fine.symbol(PI = PI, LL = LL, PL = PL)
      if(substr(fines,1,1) == "M"){
        symbol <- "SM"
      } else if(fines == "CL" || fines == "CH"){
        symbol <- "SC"
      } else if(fines == "CL-ML"){
        symbol <- "SC-SM"
      }
    }
  }

  return(symbol)
}




##------------------------------------------------------------------------------------------
##  2c. USCS SOIL CLASSIFICATION
##
##  Description:  This is a master function that calculates a soil's classification
##                according to the Unified Soil Classification System (USCS).
##                The functions A.line, USCS.coarse.symbol, and USCS.fine.symbol
##                are called by this function.
##
##  Output:   A two-element list providing a soil's USCS group symbol and name
##            symbol = USCS group symbol
##            name = USCS group name
##
##  Inputs:   pg = percent gravel (%)
##            ps = percent sand (%)
##            pf = percent fines (%)
##            Cu = coefficient of uniformity
##            Cc = coefficient of curvature
##            LL = liquid limit
##            PL = plastic limit
##            PI = plasticity index
##            sieve = vector of sieve numbers according to ASTM D422 for grain-size distribution
##            size = vector of grain sizes (in or mm) in distribution
##            percent = vector of percent passing (grain-size distribution)
##            metric = logical variable: TRUE for metric units (mm), FALSE for English units (in)
##                     [only needed if grain sizes are provided]
##
##  Notes: o  Data on the soil's grain-size distribution are required if
##            the percent fines is less than or equal to 12%.  The user has three options
##            for input to this function:
##            1.  Sieve numbers (sieve); and percent passing
##            2.  Grain sizes (size); and percent passing
##            3.  Percent components (pg, ps, pf) and coefficients of uniformity
##                and curvature (Cc and Cu)
##         o  Data on the soil's fines [either (a) LL and PL or (b) LL and PI] are
##            required if the percent fines is greater than or equal to 5%.
##         o  If percent components are specified, then Cu and Cc are needed if
##            pf < 85%.

USCS <- function(pg = NA, ps = NA, pf = NA, Cc = NA, Cu = NA,
                 LL = NA, PL = NA, PI = NA, sieve = NA, size = NA,
                 percent = NA, metric = NA){

  if(is.na(pg) == TRUE && is.na(ps) == TRUE && is.na(pf) == TRUE){
    
    ##  Obtain percent components of soil
    if(all(is.na(sieve) == FALSE)){
      size <- size.from.sieve(sieve = sieve, metric = FALSE)
      metric <- FALSE
    }
    pg <- percentComponents(size = size, percent = percent, metric = metric)$pg
    ps <- percentComponents(size = size, percent = percent, metric = metric)$ps
    pf <- percentComponents(size = size, percent = percent, metric = metric)$pf
    
    ##  Obtain coefficients of uniformity and curvature
    Cu <- grainSize.coefs(size = size, percent = percent)$Cu
    Cc <- grainSize.coefs(size = size, percent = percent)$Cc

  }
    
  ##  ---------------------------------------------------------
  ##  FINE-GRAINED SOILS
  if(pf >= 50){
    
    ##  Group symbol
    symbol <- USCS.fine.symbol(PI = PI, LL = LL, PL = PL)

    ##  Low- to medium-plasticity soils
    if(LL < 50){
      if(symbol == "CL"){      ## Classify CL Type soils		
        if(pf >= 70){
          if(pf >= 85){
            name <- "Lean Clay"
          } else if(pf  >= 70 && pf <= 85){
            if(ps >= pg){
              name <- "Lean clay with sand"
            } else if(ps < pg){
              name <- "Lean clay with gravel"
            }
          }
        } else if(pf >= 50 && pf < 70){
          if(ps >= pg){
            if(pg < 15){
              name <- "Sandy lean clay"
            } else if(pg >= 15){
              name <- "Sandy lean clay with gravel"
            }
            
          } else if(ps < pg){
            if(ps < 15){
              name <- "Gravelly lean clay"
            } else if(ps >= 15){
              name <- "Gravelly lean clay with sand"
            }
          }
        }
        
      } else if(symbol == "CL-ML"){        ## Classify CL-ML Type soils
        if(pf >= 70){
          if(pf >= 85){
            name <- "Silty clay"
          } else if(pf >= 70 && pf < 85){
            if(ps >= pg){
              name <- "Silty clay with sand"
            } else  if(ps < pg){
              name <- "Silty clay with gravel"
            }
          }
        } else if(pf >= 50 && pf < 70){
          if(ps >= pg){
            if(pg < 15){
              name <- "Sandy silty clay"
            } else if(pg >= 15){
              name <- "Sandy silty clay with gravel"
            }
          } else if(ps < pg){
            if(ps < 15){
              name <- "Gravelly silty clay"
            } else if(ps >= 15){
              name <- "Gravelly silty clay with sand"
            }
          }
        }
        
      } else if(symbol == "ML"){        ## Classify ML Type soils
        if(pf >= 70){
          if(pf >= 85){
            name <- "Silt"
          } else if(pf >= 70 && pf < 85){
            if(ps >= pg){
              name <- "Silt with sand"
            } else if(ps < pg){
              name <- "Silt with gravel"
            }
          }
        } else if(pf >= 50 && pf < 70){
          if(ps >= pg){
            if(pg < 15){
              name <- "Sandy silty"
            } else if(pg >= 15){
              name <- "Sandy silty with gravel"
            }
          } else if(ps<pg){
            if(ps < 15){
              name <- "Gravelly silt"
            } else if(ps >= 15){
              name <- "Gravelly silt with sand"
            }
          }
        }
      }

      ## High-plasticity soils
    } else if(LL >= 50){
      if(symbol == "CH"){       ## Classify CH Type soils
        if(pf >= 70){
          if(pf >= 85){
            name <- "Fat clay"
          } else if(pf >= 70 && pf < 85){
            if(ps >= pg){
              name <- "Fat clay with sand"
            } else if(ps<pg){
              name <- "Fat clay with gravel"
            }
          }
          
        } else if(pf >= 50 && pf < 70){
          if(ps >= pg){
            if(pg < 15){
              name <- "Sandy fat clay"
            }else if(pg >= 15){
              name <- "Sandy fat clay with gravel"
            }
          } else if(ps < pg){
            if(ps < 15){
              name <- "Gravelly fat clay"
            } else if(ps >= 15){
              name <- "Gravelly fat clay with sand"
            }
          }
        }
        
      } else if(symbol == "MH"){        ## Classify MH Type soils
        if(pf >= 70){
          if(pf >= 85){
            name <- "Elastic silt"
          } else if(pf >= 70 && pf < 85){
            if(ps >= pg){
              name <- "Elastic silt with sand"
            } else if(ps < pg){
              name <- "Elastic silt with gravel"
            }
          }
        } else if(pf >= 50 && pf < 70){
          if(ps >= pg){
            if(pg < 15){
              name <- "Sandy elastic silt"
            } else if(pg >= 15){
              name <- "Sandy elastic silt with gravel"
            }
          } else if(ps < pg){
            if(ps < 15){
              name <- "Gravelly elastic silt"
            } else if(ps >= 15){
              name <- "Gravelly elastic silt with sand"
            }
          }
        }
      }
    }
    
    ##  ---------------------------------------------------------
    ##  COARSE-GRAINED SOILS
  } else if(pf < 50){
    

    symbol <- USCS.coarse.symbol(pg = pg, ps = ps, pf = pf,
                                 Cc = Cc, Cu = Cu, PI = PI, LL = LL,
                                 PL = PL)
    
    ##  Group name for gravel soils
    if(symbol == "GW"){
      if(ps < 15){
        name <- "Well-graded gravel"
      } else if(ps >= 15){
        name <- "Well-graded gravel with sand"
      }
    } else if(symbol == "GP"){
      if(ps < 15){
        name <- "Poorly-graded gravel"
      } else if(ps >= 15){
        name <- "Poorly-graded gravel with sand"
      }
    } else if(symbol == "GW-GM"){
      if(ps < 15){
        name <- "Well-graded gravel with silt"
      } else if(ps >= 15){
        name <- "Well-graded gravel with silt and sand"
      }
    } else if(symbol == "GW-GC"){
      if(ps < 15){
        name <- "Well-graded gravel with clay (or silty clay)"
      } else if(ps >= 15){
        name <- "Well-graded gravel with clay & sand (or silty clay & sand)"
      }
    } else if(symbol == "GP-GM"){
      if(ps < 15){
        name <- "Poorly-graded gravel with silt"
      } else if(ps >= 15){
        name <- "Poorly-graded gravel with silt and sand"
      }
    } else if(symbol == "GP-GC"){
      if(ps < 15){
        name <- "Poorly-graded gravel with clay (or silty clay)"
      } else if(ps >= 15){
        name <- "Poorly-graded gravel with clay & sand (or silty clay & sand)"
      }
    } else if(symbol == "GM"){
      if(ps < 15){
        name <- "Silty gravel"
      } else if(ps >= 15){
        name <- "Silty gravel with sand"
      }
    } else if(symbol == "GC"){
      if(ps < 15){
        name <- "Clayey gravel"
      } else if(ps >= 15){
        name <- "Clayey gravel with sand"
      }
    } else if(symbol == "GC-GM"){
      if(ps < 15){
        name <- "Silty, clayey gravel"
      } else if(ps >= 15){
        name <- "Silty, clayey gravel with sand"
      }
      
      ##  Group name for sand soils
    } else if(symbol == "SW"){
      if(pg < 15){
        name <- "Well-graded sand"
      } else if(pg >= 15){
        name <- "Well-graded sand with gravel"
      }
    } else if(symbol == "SP"){
      if(pg < 15){
        name <- "Poorly-graded sand"
      } else if(pg >= 15){
        name <- "Poorly-graded sand with gravel"
      }
    } else if(symbol == "SW-SM"){
      if(pg < 15){
        name <- "Well-graded sand with silt"
      } else if(pg >= 15){
        name <- "Well-graded sand with silt and gravel"
      }
    } else if(symbol == "SW-SC"){
      if(pg < 15){
        name <- "Well-graded sand with clay (or silty clay)"
      } else if(pg >= 15){
        name <- "Well-graded sand with clay & gravel (or silty clay & gravel)"
      }
    } else if (symbol == "SP-SM"){
      if(pg < 15){
        name <- "Poorly-graded sand with silt"
      } else if(pg >= 15){
        name <- "Poorly-graded sand with silt and gravel"
      }
    } else if(symbol == "SP-SC"){
      if(pg < 15){
        name <- "Poorly-graded sand with clay (or silty clay)"
      } else if(pg >= 15){
        name <- "Poorly-graded sand with clay & gravel (or silty clay & gravel)"
      }
    } else if(symbol == "SM"){
      if(pg < 15){
        name <- "Silty sand"
      } else if(pg >= 15){
        name <- "Silty sand with gravel"
      }
    } else if(symbol == "SC"){
      if(pg < 15){
        name <- "Clayey sand"
      } else if(pg >= 15){
        name <- "Clayey sand with gravel"
      }
    } else if(symbol == "SC-SM"){
      if(pg < 15){
        name <- "Silty, clayey sand"
      } else if(pg >= 15){
        name <- "Silty, clayey sand with gravel"
      }
    }
  }
  ## Return soil group symbol and name
  return(list(symbol = symbol, name = name))
}

##  Example
##  USCS(pg = 15, ps = 34, pf = 51, Cc = 1, Cu = 4, LL = 40, PL = 10)
