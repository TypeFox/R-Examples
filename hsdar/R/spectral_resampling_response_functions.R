get_response_function <- function(sensor)
{
  pc <- if (nchar(deparse(sys.calls()[[sys.nframe()-1]])) > 22) 
    substr(deparse(sys.calls()[[sys.nframe()-1]]),1, 22) != "list.available.sensors" else TRUE
  response <- switch(sensor,
                     "RapidEye"=get_RapidEye_response(),
                     "WorldView2-8"=get_wv2_8_response(pc),
                     "Quickbird"=get_quickbird_response(pc),                 
                     "WorldView2-4"=get_wv2_4_response(pc),
                     "Landsat4"=get_landsat4_response(),
                     "Landsat5"=get_landsat5_response(),
                     "Landsat7"=get_landsat7_response(),
                     "Landsat8"=get_landsat8_response(),
#                      "Modis"=get_TerraModis_response(),
                     NULL
              )
  return(response)
}

get_RapidEye_response <- function ()
{
  RapidEye_response <- NULL
  rm(RapidEye_response)
  data("RapidEye_response", package = "hsdar", envir = environment())
  response <- RapidEye_response
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 419
  attr(response, "maxwl")    <- 901
  attr(response, "stepsize") <- 1
  return(response)
}

get_wv2_8_response <- function (pc)
{
  WV_2_8_response <- NULL
  rm(WV_2_8_response)
  
  ## Copyright by DigitalGlobe, Inc. All Rights Reserved
  data("WV_2_8_response", package = "hsdar", envir = environment())
  
  response <- WV_2_8_response
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 349
  attr(response, "maxwl")    <- 1101
  attr(response, "stepsize") <- 1
  if (pc)
    cat("Copyright of spectral response function by DigitalGlobe, Inc. All Rights Reserved\n")
  return(response)
}

get_wv2_4_response <- function (pc)
{
  WV_2_8_response <- NULL
  rm(WV_2_8_response)
  
  ## Copyright by DigitalGlobe, Inc. All Rights Reserved
  data("WV_2_8_response", package = "hsdar", envir = environment())
  
  response <- WV_2_8_response[,c(2,4,5,7)]
  names(response)[4] <- "NIR"
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 349
  attr(response, "maxwl")    <- 1101
  attr(response, "stepsize") <- 1
  if (pc)
    cat("Copyright of spectral response function by DigitalGlobe, Inc. All Rights Reserved\n")
  return(response)
}

get_quickbird_response <- function (pc)
{
  Quickbird_response <- NULL
  rm(Quickbird_response)
  
  ## Copyright by DigitalGlobe, Inc. All Rights Reserved
  data("Quickbird_response", package = "hsdar", envir = environment())
  
  response <- Quickbird_response
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 300
  attr(response, "maxwl")    <- 1100
  attr(response, "stepsize") <- 5
  if (pc)
    cat("Copyright of spectral response function by DigitalGlobe, Inc. All Rights Reserved\n")
  return(response)
}

get_landsat4_response <- function ()
{
  Landsat_4_response <- NULL
  rm(Landsat_4_response)
  data("Landsat_4_response", package = "hsdar", envir = environment())
  response <- Landsat_4_response
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 410
  attr(response, "maxwl")    <- 2400
  attr(response, "stepsize") <- 1
  return(response)
}

get_landsat5_response <- function ()
{
  Landsat_5_response <- NULL
  rm(Landsat_5_response)
  data("Landsat_5_response", package = "hsdar", envir = environment())
  response <- Landsat_5_response
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 418
  attr(response, "maxwl")    <- 2401
  attr(response, "stepsize") <- 1
  return(response)
}

get_landsat7_response <- function ()
{
  Landsat_7_response <- NULL
  rm(Landsat_7_response)
  data("Landsat_7_response", package = "hsdar", envir = environment())
  response <- Landsat_7_response
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 434
  attr(response, "maxwl")    <- 2401
  attr(response, "stepsize") <- 1
  return(response)
}

get_landsat8_response <- function ()
{
  Landsat_8_response <- NULL
  rm(Landsat_8_response)
  data("Landsat_8_response", package = "hsdar", envir = environment())
  response <- Landsat_8_response
  attr(response, "wlunit")   <- "nm"
  attr(response, "minwl")    <- 427
  attr(response, "maxwl")    <- 2355
  attr(response, "stepsize") <- 1
  return(response)
}

# 
# get_TerraModis_response <- function ()
# {
#   Landsat_7_response <- NULL
#   rm(Landsat_7_response)
#   data("Landsat_7_response", package = "hsdar")
#   response <- Landsat_7_response
#   attr(response, "wlunit")   <- "nm"
#   attr(response, "minwl")    <- 434
#   attr(response, "maxwl")    <- 2401
#   attr(response, "stepsize") <- 1
#   return(response)
# }