# This file is a part of the helsinki package (http://github.com/rOpenGov/helsinki)
# in association with the rOpenGov project (ropengov.github.io)

# Copyright (C) 2010-2014 Juuso Parkkinen, Leo Lahti and Joona Lehtomaki / Louhos <louhos.github.com>. 
# All rights reserved.

# This program is open source software; you can redistribute it and/or modify 
# it under the terms of the FreeBSD License (keep this notice): 
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Retrieve data from Helsinki Region Environmental Services
#'
#' Retrieves data from Helsinki Region Environmental
#' Services Authority (Helsingin seudun ymparistopalvelu HSY) 
#' http://www.hsy.fi/seututieto/kaupunki/paikkatiedot/Sivut/Avoindata.aspx
#' For data description (in Finnish) see:
#' http://www.hsy.fi/seututieto/Documents/Paikkatiedot/Tietokuvaukset_kaikki.pdf. 
#' The data copyright (C) HSY 2011.
#'
#' @param which.data A string. Specify the name of the retrieved HSY data set. Options: Vaestotietoruudukko; Rakennustietoruudukko; SeutuRAMAVA_kosa; SeutuRAMAVA_tila. These are documented in HSY data description document (see above).
#' @param which.year An integer. Specify the year for the data to be retrieved.
#' @param data.dir A string. Specify a temporary folder for storing downloaded data.
#' @param verbose logical. Should R report extra information on progress? 
#'
#' @return Shape object (from SpatialPolygonsDataFrame class)
#' @import maptools
#' @importFrom sp CRS
#' @export
#' @references See citation("helsinki") 
#' @author Juuso Parkkinen and Leo Lahti \email{louhos@@googlegroups.com}
#' @examples vaesto.sp <- get_hsy("Vaestotietoruudukko")
#' @keywords utilities

get_hsy <- function (which.data=NULL, which.year=2013, data.dir=tempdir(), verbose=TRUE) {
  
  message("IMPORTANT NOTE! HSY open data services have been recently updated and get_hsy() function is outdated! It will be updated soon, meanwhile use the services directly at https://www.hsy.fi/fi/asiantuntijalle/avoindata/Sivut/default.aspx.")
  return(NULL)
  
  if (is.null(which.data)) {
    message("Available HSY datasets:
  'Vaestotietoruudukko': Ruutukohtaista tietoa vaeston lukumaarasta, ikajakaumasta ja asumisvaljyydesta. Vuodet: 1997-2003, 2008-2013.
  'Rakennustietoruudukko': Ruutukohtaista tietoa rakennusten lukumaarasta, kerrosalasta, kayttotarkoituksesta ja aluetehokkuudesta. Ruutukoko 500x500 metria. Vuodet: 1997-2003, 2008-2013.
  'SeutuRAMAVA_kosa': kaupunginosittain summattua tietoa rakennusmaavarannosta. Vuodet: 2010-2013.
  'SeutuRAMAVA_tila': tilastoalueittain summattua tietoa rakennusmaavarannosta. Vuodet: 1997-2001, 2010-2013")
    stop("Please specify 'which.data'")
  }
  
  if (which.data=="SeutuRAMAVA_tila" & which.year==2012)
    stop("Data for SeutuRAMAVA_tila 2012 is not readable with maptools::readShapePoly. You can try using the rgdal-package if you need to access that data!")
  
  # Create data.dir if it does not exist
  if (!file.exists(data.dir))
    dir.create(data.dir)
  
  
  ## Download data ----------------------------------------------------
  
  # Specify download url
  if (which.data=="Vaestotietoruudukko") {
    if (which.year==2013)
      zip.file <- paste0("Vaestotietoruudukko", "_", which.year, "_SHP.zip")
    else if (which.year==2012)
      zip.file <- paste0("Vaestoruudukko", "_", which.year, ".zip")
    else
      zip.file <- paste0("Vaestoruudukko", "_SHP_", which.year, ".zip")
    
  } else if (which.data=="Rakennustietoruudukko") {
    if (which.year==2013)
      zip.file <- paste0("Rakennustietoruudukko", "_", which.year, "_SHP.zip")
    else if (which.year==2012)
      zip.file <- paste0("Rakennustietoruudukko", "_", which.year, ".zip")
    else
      zip.file <- paste0("Rakennustietoruudukko", "_SHP_", which.year, ".zip")
    
  } else if (which.data=="SeutuRAMAVA_kosa") {
    if (which.year==2013) {
      zip.file <- paste0("SeutuRamava_kosa_", which.year, "_SHP.zip")
    } else if (which.year %in% c(2011, 2012)) {
      zip.file <- paste0("SeutuRAMAVA_SHP_", which.year, ".zip")
    } else {
      zip.file <- paste0("SeutuRAMAVA_SHP.zip")
    }
    
  } else if (which.data=="SeutuRAMAVA_tila") {
    if (which.year==2013) {
      zip.file <- paste0("SeutuRamava_tila_", which.year, "_SHP.zip")
    } else if (which.year==2012) {
      zip.file <- paste0("SeutuRAMAVA_tila_shp", which.year, ".zip")
    } else {
      zip.file <- paste0("SeutuRAMAVA_tila_", which.year, "_SHP.zip")
    }
  } else {
    stop("Invalid 'which.data' argument")
  }
  
  # Download data
  if (which.data=="SeutuRAMAVA_kosa" & which.year==2010) {
    remote.zip <- paste0("http://www.hsy.fi/seututieto/Documents/Paikkatiedot/", zip.file) 
  } else {
    remote.zip <- paste0("http://www.hsy.fi/seututieto/kaupunki/paikkatiedot/Documents/", zip.file)
  }
  local.zip <- file.path(data.dir, zip.file)
  if (!file.exists(local.zip)) {
    # Check whether url available
    if (!RCurl::url.exists(remote.zip)) {
      message(paste("Sorry! Url", remote.zip, "not available!\nReturned NULL."))
      return(NULL)
    }
    
    if (verbose)
      message("Dowloading ", remote.zip, "\ninto ", local.zip, "\n")
    utils::download.file(remote.zip, destfile = local.zip, quiet=!verbose)
  } else {
    if (verbose)
      message("File ", local.zip, " already found, will not download again!")
  }  
  
  ## Process data -----------------------------------------------
  
  # Unzip the downloaded zip file
  utils::unzip(local.zip, exdir = data.dir)
  
  # Define shapefile
  if (which.data=="Vaestotietoruudukko") {
    if (verbose)
      message("For detailed description of ", which.data, " see\nhttp://www.hsy.fi/seututieto/kaupunki/paikkatiedot/Documents/Vaestoruudukko.pdf")
    if (which.year %in% c(2012, 2013))
      sp.file <- paste0(data.dir, "/Vaestotietoruudukko_", which.year, ".shp")
    else
      sp.file <- paste0(data.dir, "/Vaestoruudukko_", which.year, ".shp")
    
  } else if (which.data=="Rakennustietoruudukko") {
    if (verbose)
      message("For detailed description of ", which.data, " see\nhttp://www.hsy.fi/seututieto/kaupunki/paikkatiedot/Documents/Rakennustietoruudukko.pdf")
    if (which.year==2012)
      sp.file <- paste0(data.dir, "/Rakennustitetoruudukko_", which.year, ".shp")
    else
      sp.file <- paste0(data.dir, "/Rakennustietoruudukko_", which.year, ".shp")
    
  } else if (which.data=="SeutuRAMAVA_kosa") {
    if (verbose) 
      message("For detailed description of ", which.data, " see\nhttp://www.hsy.fi/seututieto/kaupunki/paikkatiedot/Documents/SeutuRAMAVA.pdf")
    if (which.year==2013)
      sp.file <- paste0(data.dir, "/SeutuRamava_kosa_", which.year, ".shp")
    else if (which.year==2012)
      sp.file <- paste0(data.dir, "/SeutuRAMAVA_kosa_", which.year, ".shp")
    else
      sp.file <- paste0(data.dir, "/SeutuRAMAVA_", which.year, ".shp")
    
  } else if (which.data=="SeutuRAMAVA_tila") {
    if (verbose)
      message("For detailed description of ", which.data, " see\nhttp://www.hsy.fi/seututieto/kaupunki/paikkatiedot/Documents/SeutuRAMAVA.pdf")
    if (which.year==2013)
      sp.file <- paste0(data.dir, "/SeutuRamava_tila_", which.year, ".shp")
    else if (which.year==2012)
      sp.file <- paste0(data.dir, "/SeutuRAMAVA_Tila_", which.year, ".shp")
    else
      sp.file <- paste0(data.dir, "/SeutuRAMAVA_", which.year, "_SHP.shp")
  }
  
  # Read shapefile and add coordinate information manually (ETRS-GK25 -> EPSG:3879)
  p4s <- "+init=epsg:3879 +proj=tmerc +lat_0=0 +lon_0=25 +k=1 +x_0=25500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  sp <- maptools::readShapePoly(fn=sp.file, proj4string=sp::CRS(p4s))
    
  # Add KATAKER to rakennustieto, mailaa '11' SePe:lle
  if (which.data=="Rakennustietoruudukko") {
    KATAKER.key <- kataker_key()
    kk.df <- data.frame(list(KATAKER = as.integer(names(KATAKER.key)), description = KATAKER.key))
    temp <- sp@data
    temp$KATAKER1.description <- kk.df$description[match(temp$KATAKER1, kk.df$KATAKER)]
    temp$KATAKER2.description <- kk.df$description[match(temp$KATAKER2, kk.df$KATAKER)]
    temp$KATAKER3.description <- kk.df$description[match(temp$KATAKER3, kk.df$KATAKER)]
    sp@data <- temp
    message("\nAdded KATAKER descriptions to Rakennustietoruudukko")
  }
  
  # Fix characters in SeutuRAMAVA data
  if (which.data=="SeutuRAMAVA_kosa") {
    for (nam in c("OMLAJI_1S", "OMLAJI_2S", "OMLAJI_3S", "NIMI", "NIMI_SE"))    
      sp[[nam]] <-  factor(iconv(sp[[nam]], from = "latin1", to = "UTF-8"))
  }
  if (which.data=="SeutuRAMAVA_tila") {
    for (nam in c("OMLAJI_1S", "OMLAJI_2S", "OMLAJI_3S", "NIMI"))    
      sp[[nam]] <-  factor(iconv(sp[[nam]], from = "latin1", to = "UTF-8"))
  }
  
  if (verbose)
    message("\nData loaded succesfully!")
  return(sp)
}


kataker_key <- function () {
  KATAKER.key <- c(
    "11"   = "Yhden asunnon talot", 
    "12"  = "Kahden asunnon talot", 
    "13"  = "Muut erilliset pientalot", 
    "21"  = "Rivitalot", 
    "22"  = "Ketjutalot", 
    "32"  = "Luhtitalot", 
    "39"  = "Muut kerrostalot", 
    "41"  = "Vapaa-ajan asunnot",
    "111" = "Myymalahallit", 
    "112" = "Liike- ja tavaratalot, kauppakeskukset",
    "119" = "Myymalarakennukset ", 
    "121" = "Hotellit, motellit, matkustajakodit, kylpylahotellit", 
    "123" = "Loma- lepo- ja virkistyskodit", 
    "124" = "Vuokrattavat lomamokit ja osakkeet (liiketoiminnallisesti)", 
    "129" = "Muut majoitusliikerakennukset",
    "131" = "Asuntolat, vanhusten palvelutalot, asuntolahotellit",
    "139" = "Muut majoitusrakennukset", 
    "141" = "Ravintolat, ruokalat ja baarit", 
    "151" = "Toimistorakennukset", 
    "161" = "Rautatie- ja linja- autoasemat, lento- ja satamaterminaalit", 
    "162" = "Kulkuneuvojen suoja- ja huoltorakennukset", 
    "163" = "Pysakointitalot", 
    "164" = "Tietoliikenteen rakennukset", 
    "165" = "Muut liikenteen rakennukset", 
    "169" = "Muut liikenteen rakennukset", 
    "211" = "Keskussairaalat", 
    "213" = "Muut sairaalat", 
    "214" = "Terveyskeskukset", 
    "215" = "Terveydenhoidon erityislaitokset (mm. kuntoutuslaitokset)", 
    "219" = "Muut terveydenhoitorakennukset", 
    "221" = "Vanhainkodit", 
    "222" = "Lastenkodit, koulukodit", 
    "223" = "Kehitysvammaisten hoitolaitokset", 
    "229" = "Muut huoltolaitosrakennukset", 
    "231" = "Lasten paivakodit", 
    "239" = "Muut sosiaalitoimen rakennukset", 
    "241" = "Vankilat", 
    "311" = "Teatterit, konsertti- ja kongressitalot, oopperat", 
    "312" = "Elokuvateatterit",
    "322" = "Kirjastot", 
    "323" = "Museot, taidegalleriat",
    "324" = "Nayttelyhallit", 
    "331" = "Seurain-, nuoriso- yms. talot",
    "341" = "Kirkot, kappelit, luostarit, rukoushuoneet",
    "342" = "Seurakuntatalot", 
    "349" = "Muut uskonnollisten yhteisojen rakennukset", 
    "351" = "Jaahallit", 
    "352" = "Uimahallit", 
    "353" = "Tennis-, squash- ja sulkapallohallit",
    "354" = "Monitoimi- ja muut urheiluhallit",
    "359" = "Muut urheilu- ja kuntoilurakennukset", 
    "369" = "Muut kokoontumis- rakennukset", 
    "511" = "Peruskoulut, lukiot ja muut", 
    "521" = "Ammatilliset oppilaitokset", 
    "531" = "Korkeakoulu- rakennukset",
    "532" = "Tutkimuslaitosrakennukset", 
    "541" = "Jarjestojen, liittojen, tyonantajien yms.  opetusrakennukset", 
    "549" = "Muualla luokittelemattomat opetusrakennukset", 
    "611" = "Voimalaitosrakennukset", 
    "613" = "Yhdyskuntatekniikan rakennukset", 
    "691" = "Teollisuushallit", 
    "692" = "Teollisuus- ja pienteollisuustalot", 
    "699" = "Muut teollisuuden tuotantorakennukset", 
    "711" = "Teollisuusvarastot",
    "712" = "Kauppavarastot", 
    "719" = "Muut varastorakennukset",
    "721" = "Paloasemat", 
    "722" = "Vaestonsuojat", 
    "723" = "Halytyskeskukset",
    "729" = "Muut palo- ja pelastustoimen rakennukset", 
    "811" = "Navetat, sikalat, kanalat yms.", 
    "819" = "Elainsuojat, ravihevostallit, maneesit", 
    "891" = "Viljankuivaamot ja viljan sailytysrakennukset, siilot", 
    "892" = "Kasvihuoneet", 
    "893" = "Turkistarhat", 
    "899" = "Muut maa-, metsa- ja kalatalouden rakennukset", 
    "931" = "Saunarakennukset",
    "941" = "Talousrakennukset", 
    "999" = "Muut rakennukset", 
    "999999999" = "Puuttuvan tiedon merkki")
  
  KATAKER.key
}
