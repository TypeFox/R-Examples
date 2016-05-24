#' Create a DRAC input data template (v1.1)
#'
#' This function returns a DRAC input template (v1.1) to be used in conjunction
#' with the use_DRAC() function
#' 
#' @param nrow \code{\link{integer}} (with default): specifies the number of rows
#' of the template (i.e., the number of data sets you want to submit)
#' 
#' @param notification \code{\link{logical}} (with default): show or hide the
#' notification
#'
#' @return A list.
#' 
#' @author Christoph Burow, University of Cologne (Germany)
#'
#' @references
#'
#' Durcan, J.A., King, G.E., Duller, G.A.T., 2015. DRAC: Dose Rate and Age Calculator for trapped charge dating.
#' Quaternary Geochronology 28, 54-61. doi:10.1016/j.quageo.2015.03.012
#'
#' @seealso \code{\link{as.data.frame}} \code{\link{list}} 
#'
#' @examples
#' 
#' # create a new DRAC input input
#' input <- template_DRAC()
#' 
#' # show content of the input
#' print(input)
#' print(input$`Project ID`)
#' print(input[[4]])
#' 
#' 
#' ## Example: DRAC Quartz example
#' # note that you only have to assign new values where they 
#' # are different to the default values
#' input$`Project ID` <- "DRAC-Example"
#' input$`Sample ID` <- "Quartz"
#' input$`Conversion factors` <- "AdamiecAitken1998"
#' input$`ExternalU (ppm)` <- 3.4
#' input$`errExternal U (ppm)` <- 0.51
#' input$`External Th (ppm)` <- 14.47
#' input$`errExternal Th (ppm)` <- 1.69
#' input$`External K (%)` <- 1.2
#' input$`errExternal K (%)` <- 0.14
#' input$`Calculate external Rb from K conc?` <- "N"
#' input$`Calculate internal Rb from K conc?` <- "N"
#' input$`Scale gammadoserate at shallow depths?` <- "N"
#' input$`Grain size min (microns)` <- 90
#' input$`Grain size max (microns)` <- 125
#' input$`Water content ((wet weight - dry weight)/dry weight) %` <- 5
#' input$`errWater content %` <- 2
#' input$`Depth (m)` <- 2.2
#' input$`errDepth (m)` <- 0.22
#' input$`Overburden density (g cm-3)` <- 1.8
#' input$`errOverburden density (g cm-3)` <- 0.1
#' input$`Latitude (decimal degrees)` <- 30.0000
#' input$`Longitude (decimal degrees)` <- 70.0000
#' input$`Altitude (m)` <- 150
#' input$`De (Gy)` <- 20
#' input$`errDe (Gy)` <- 0.2
#' 
#' # use DRAC
#' \dontrun{
#' output <- use_DRAC(input)
#' }
#' 
#' @export
template_DRAC <- function(nrow = 1, notification = TRUE) {
  
  ## TODO:
  # 1 - allow mineral specific presets; new argument 'mineral'
  # 2 - add option to return the DRAC example data set
  
  if (nrow < 0 | nrow > 33) 
    stop("'nrow' must be a number between 0 and 33.", call. = FALSE)
  
  ## LEGAL NOTICE ----
  messages <- list("\n",
                   "\t-------------------- IMPORTANT NOTE ------------------------\n",
                   "\t This function returns a DRAC input template to be used in ",
                   "\t conjunction with the use_DRAC() function.  \n",
                   "\t The template was reproduced with great care, but we do not",
                   "\t take any responsibility and we are not liable for any ",
                   "\t mistakes or unforeseen misbehaviour.",
                   "\t Note that this template is only compatible with DRAC",
                   "\t version 1.1. Before using this template make sure that",
                   "\t this is the correct version, otherwise expect unspecified",
                   "\t errors.\n",
                   "\t Please ensure you cite the use of DRAC in your work,",
                   "\t published or otherwise. Please cite the website name and",
                   "\t version (e.g. DRAC v1.1) and the accompanying journal",
                   "\t article:",
                   "\t Durcan, J.A., King, G.E., Duller, G.A.T., 2015.",
                   "\t DRAC: Dose rate and age calculation for trapped charge",
                   "\t dating. Quaternary Geochronology 28, 54-61. \n",
                   "\t Set 'notification = FALSE' to hide this message. \n",
                   "\t-------------------- IMPORTANT NOTE ------------------------",
                   "\n")
  
  if (notification) lapply(messages, message)
  
  # CREATE TEMPLATE ----
  template <- list(
    
    `Project ID` = 
      structure(rep("RLum", nrow), required = TRUE, allowsX = FALSE, key = "TI:1",
                description = "Inputs can be alphabetic, numeric or selected symbols (/ - () [] _). Spaces are not permitted."), # 
    
    `Sample ID` = 
      structure(rep("999", nrow), required = TRUE, allowsX = FALSE, key = "TI:2",
                description = "Inputs can be alphabetic, numeric or selected symbols (/ - () [] _). Spaces are not permitted."), #
    
    `Mineral` = 
      structure(factor(rep("Q", nrow), c("Q", "F", "PM")), required = TRUE, allowsX = FALSE, key = "TI:3",
                description = "The mineral used for dating: quartz, feldspar or polymineral. Input must be 'Q', 'F' or 'PM'."), #
    
    `Conversion factors` = 
      structure(factor(rep("Liritzisetal2013", nrow), c("AdamiecAitken1998", "Guerinetal2011", "Liritzisetal2013", "X")), required = FALSE, allowsX = TRUE, key = "TI:4",
                description = "The conversion factors required to calculate dose rates from radionuclide concentrations. Users have the option of datasets from Adamiec and Aitken (1998), Guerin et al. (2011) or Liritzis et al. (2013). Input must be 'AdamiecAitken1998', 'Guerinetal2011', 'Liritzisetal2013' or 'X' if conversion factors are not required."), #
    
    `ExternalU (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:5",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), # 
    
    `errExternal U (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:6",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `External Th (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:7",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `errExternal Th (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:8",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `External K (%)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:9",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `errExternal K (%)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:10",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `External Rb (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:11",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `errExternal Rb (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:12",
                description = "Radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `Calculate external Rb from K conc?` = 
      structure(factor(rep("Y", nrow), c("Y", "N")), required = FALSE, allowsX = FALSE, key = "TI:13",
                description = "Option to calculate a Rubidium concentration from Potassium, using the 270:1 ratio suggested by Mejdahl (1987). Input should be yes 'Y' or no 'N'."), #
    
    `Internal U (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:14",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `errInternal U (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:15",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `Internal Th (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:16",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `errInternal Th (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:17",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `Internal K (%)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:18",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `errInternal K (%)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:19",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `Rb (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:20",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `errRb (ppm)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:21",
                description = "Internal radionuclide concentrations in parts per million for Uranium, Thorium and Rubidium and % for Potassium. Inputs must be 0 or positive and should not be left blank."), #
    
    `Calculate internal Rb from K conc?` = 
      structure(factor(rep("Y", nrow), c("Y", "N", "X")), required = FALSE, allowsX = TRUE, key = "TI:22",
                description = "Option to calculate an internal Rubidium concentration from Potassium, using the 270:1 ratio suggested by Mejdahl (1987). Input should be yes 'Y' or no 'N'."), #
    
    `User external alphadoserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:23",
                description = "Users may input directly measured values for external alpha, beta and gamma dose rates (in Gy.ka-1). Any positive inputs in these fields will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and should not be left blank"), #
    
    `errUser external alphadoserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:24",
                description = "Users may input directly measured values for external alpha, beta and gamma dose rates (in Gy.ka-1). Any positive inputs in these fields will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and should not be left blank"), #
    
    `User external betadoserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:25",
                description = "Users may input directly measured values for external alpha, beta and gamma dose rates (in Gy.ka-1). Any positive inputs in these fields will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and should not be left blank"), #
    
    `errUser external betadoserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:26",
                description = "Users may input directly measured values for external alpha, beta and gamma dose rates (in Gy.ka-1). Any positive inputs in these fields will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and should not be left blank"), #
    
    `User external gamma doserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:27",
                description = "Users may input directly measured values for external alpha, beta and gamma dose rates (in Gy.ka-1). Any positive inputs in these fields will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and should not be left blank"), #
    
    `errUser external gammadoserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:28",
                description = "Users may input directly measured values for external alpha, beta and gamma dose rates (in Gy.ka-1). Any positive inputs in these fields will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and should not be left blank"), #
    
    `User internal doserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:29",
                description = "Users may input an internal dose rate (either alpha, beta or the sum of the two; in Gy.ka-1). DRAC will assume that this value has already been corrected for attenuation. Inputs in this field will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and not left blank."), #
    
    `errUser internal doserate (Gy.ka-1)` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:30",
                description = "Users may input an internal dose rate (either alpha, beta or the sum of the two; in Gy.ka-1). DRAC will assume that this value has already been corrected for attenuation. Inputs in this field will override dose rates calculated from radionuclide concentrations. Inputs should be 0 or positive and not left blank."), #
    
    `Scale gammadoserate at shallow depths?` = 
      structure(factor(rep("Y", nrow), c("Y", "N")), required = FALSE, allowsX = FALSE, key = "TI:31",
                description = "Users may choose to scale gamma dose rates for samples taken within 0.3 m of the ground surface. The scaling factors of Aitken (1985) are used. Input should be yes 'Y' or no 'N'."), #
    
    `Grain size min (microns)` = 
      structure(rep(100, nrow), required = TRUE, allowsX = FALSE, key = "TI:32",
                description = "The grain size range analysed. DRAC can be used for the grain size ranges between 1 and 1000 microns. Inputs should range between 1 and 1000 and not be left blank."), #
    
    `Grain size max (microns)` = 
      structure(rep(150, nrow), required = TRUE, allowsX = FALSE, key = "TI:33",
                description = "The grain size range analysed. DRAC can be used for the grain size ranges between 1 and 1000 microns. Inputs should range between 1 and 1000 and not be left blank."), #
    
    `alpha-Grain size attenuation` = 
      structure(factor(rep("Brennanetal1991", nrow), c("Bell1980", "Brennanetal1991")), required = TRUE, allowsX = FALSE, key = "TI:34",
                description = "The grain size attenuation factors for the alpha dose rate. Users have the option of datasets from Bell (1980) and Brennan et al. (1991). Input must be 'Bell1980' or 'Brennanetal1991'."), #
    
    `beta-Grain size attenuation ` = 
      structure(factor(rep("Guerinetal2012-Q", nrow), c("Mejdahl1979", "Brennan2003", "Guerinetal2012-Q", "Guerinetal2012-F")), required = TRUE, allowsX = FALSE, key = "TI:35",
                description = "The grain size attenuation factors for the beta dose rate. Users have the option of datasets from Mejdahl (1979), Brennan (2003) and Guerin et al. (2012) for quartz or feldspar. Input must be 'Mejdahl1979', 'Brennan2003', 'Guerinetal2012-Q' or 'Guerinetal2012-F' ."), #
    
    `Etch depth min (microns)` = 
      structure(rep(8, nrow), required = TRUE, allowsX = FALSE, key = "TI:36",
                description = "The user defined etch depth range (microns). Inputs should range between 0 and 30 and not be left blank."), #
    
    `Etch depth max (microns)` = 
      structure(rep(10, nrow), required = TRUE, allowsX = FALSE, key = "TI:37",
                description = "The user defined etch depth range (microns). Inputs should range between 0 and 30 and not be left blank."), #
    
    `beta-Etch depth attenuation factor` = 
      structure(factor(rep("Bell1979", nrow), c("Bell1979", "Brennan2003", "X")), required = FALSE, allowsX = TRUE, key = "TI:38",
                description = "The etch depth attenuation factors for the beta dose rate. Users have the option of datasets from Bell (1979) and Brennan (2003). Input must be 'Bell1979' or 'Brennan2003'. Note: only the dataset of Bell (1980) is provided for attenuation of the alpha dose rate by etching."), #
    
    `a-value` = 
      structure(rep(0, nrow), required = FALSE, allowsX = TRUE, key = "TI:39",
                description = "Alpha track efficiency value and uncertainty defined by the user. Inputs should be 0 or positive and not left blank."), #
    
    `erra-value` = 
      structure(rep(0, nrow), required = TRUE, allowsX = TRUE, key = "TI:40",
                description = "Alpha track efficiency value and uncertainty defined by the user. Inputs should be 0 or positive and not left blank."), #
    
    `Water content ((wet weight - dry weight)/dry weight) %` = 
      structure(rep(0, nrow), required = TRUE, allowsX = FALSE, key = "TI:41",
                description = "Sediment water content (%) over the burial period. Inputs should be 0 or positive and not be left blank."), #
    
    `errWater content %` = 
      structure(rep(0, nrow), required = FALSE, allowsX = FALSE, key = "TI:42",
                description = "Sediment water content (%) over the burial period. Inputs should be 0 or positive and not be left blank."), #
    
    `Depth (m)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:43",
                description = "Depth and uncertainty from which sample was extracted beneath the ground surface. Inputs should be 0 or positive and not left blank. If user defined Dc will be used then an 'X' must be input."), #
    
    `errDepth (m)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:44",
                description = "Depth and uncertainty from which sample was extracted beneath the ground surface. Inputs should be 0 or positive and not left blank. If user defined Dc will be used then an 'X' must be input."), #
    
    `Overburden density (g cm-3)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:45",
                description = "Density of the overlying sediment matrix from which the sample was taken. Inputs should be 0 or positive and not be left blank. If user defined Dc will be used then an 'X' must be input."), #
    
    `errOverburden density (g cm-3)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:46",
                description = "Density of the overlying sediment matrix from which the sample was taken. Inputs should be 0 or positive and not be left blank. If user defined Dc will be used then an 'X' must be input."), #
    
    `Latitude (decimal degrees)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:47",
                description = "Latitude and longitude of sample location (in degree decimals). Positive values should be used for northern latitudes and eastern longitudes and negative values for southern latitudes and western longitudes. Inputs should range from -90 to 90 degrees for latitudes and -180 to 180 degrees for longitude. If user defined Dc will be used then an 'X' must be input."), # 
    
    `Longitude (decimal degrees)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:48",
                description = "Latitude and longitude of sample location (in degree decimals). Positive values should be used for northern latitudes and eastern longitudes and negative values for southern latitudes and western longitudes. Inputs should range from -90 to 90 degrees for latitudes and -180 to 180 degrees for longitude. If user defined Dc will be used then an 'X' must be input."), # 
    
    `Altitude (m)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:49",
                description = "Altitude of sample location in metres above sea level. Input should be less than 5000 and not left blank. If user defined Dc will be used then an 'X' must be input."), #
    
    `User cosmicdoserate (Gy.ka-1)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:50",
                description = "Users may input a cosmic dose rate (in Gy.ka-1). Inputs in these fields will override the DRAC calculated cosmic dose rate. Inputs should be positive or 'X' if not required, and not left blank."), #
    
    `errUser cosmicdoserate (Gy.ka-1)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:51",
                description = "Users may input a cosmic dose rate (in Gy.ka-1). Inputs in these fields will override the DRAC calculated cosmic dose rate. Inputs should be positive or 'X' if not required, and not left blank."), #
    
    `De (Gy)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:52",
                description = "Sample De and uncertainty (in Gy). Inputs should be positive or 'X' if not required, and not left blank."), #
    
    `errDe (Gy)` = 
      structure(rep("X", nrow), required = FALSE, allowsX = TRUE, key = "TI:53",
                description = "Sample De and uncertainty (in Gy). Inputs should be positive or 'X' if not required, and not left blank.") #
  )   
  
  
  ## RETURN VALUE ---
  # add an additional DRAC class so we can define our own S3 method for as.data.frame
  class(template) <- c("DRAC.list", "list")
  invisible(template)
}