# see
#http://sdmx.org/wp-content/uploads/2013/09/SDMX_2_1-SECTION_07_WebServicesGuidelines_2013-04.pdf
# section 5.8 regarding errors.

require("RJSDMX")


###################### NOTE re try(..., silent=TRUE) ######################
#  silent=TRUE only works properly, using level OFF in configuration file
z <- try(getSDMX("OECD", 'G20_PRICES.CAB.CP.IXOB.M'), silent=TRUE)

##################################################

#  https://github.com/amattioc/SDMX/wiki

# to install devel version from Github
# require(devtools)
# old install_github(repo = "SDMX", username = "amattioc", subdir = "RJSDMX")
# install_github(repo = "amattioc/SDMX", subdir = "RJSDMX")
# check installed version
# installed.packages()["RJSDMX",c("Package","Version")] 
# used 1.1 testing to 6 Nov 2014
#      1.2 installed from github 6 Nov 2014
#      1.3 installed from github 11 Dec 2014
#      1.3 installed from github 22 Dec 2014
#      1.3 installed from github  5 Feb 2015
#      1.4.2 installed from github  14 Sept 2015

# Package rJava may be needed interactively for experimenting but should be found
#   in the the namespace when everything is working.
# require("rJava")


# When the package is loaded there is an indication if/where 
# a configuration file has been found.
# Its location should be set with the SDMX_CONF environment variable.
#  e.g  export SDMX_CONF=/home/paul/.SdmxClient
# Details about contents are at https://github.com/amattioc/SDMX/wiki/Configuration
# The configuration file can be used to control the level of std output about
# warnings and errors. This mostly seems to be coming directly from the java,
# rather than passed back to R (which would be more usual for R packages as it
# can then be masked in R by try(), etc, if that makes sense.)
# R users many want to set 
#SDMX.level = WARNING
#java.util.logging.ConsoleHandler.level = WARNING
# to limit output to what would more usually be expected in R sessions.


##########  Notes on finding identifiers (e.g. EuroStat) ##################

# Finding series identifiers is difficult and
#  the mneumonics are obscure. Needs documentation.

# on finding identifiers
# http://epp.eurostat.ec.europa.eu/portal/page/portal/sdmx_web_services/about_eurostat_data
# >SDMX queries tutorial
# and 
#http://epp.eurostat.ec.europa.eu/portal/page/portal/sdmx_web_services/sdmx_tutorial/ident_data

#   but real, use   sdmxHelp()


# sdmxHelp()  # This is very helpful

#>EUROSTAT  ei_nama_q : Main aggregates - quarterly
# 	eil_nama_q  > 	>FREQ: Q
#			>UNIT : MIO-EUR
#			>S_ADJ: NSA
#			>P_ADJ: CP  (current prices)
#			>INDIC: NA-B1GP  (GDP at market prices)
#			>GEO: IT  (Italy)

# This shows all ei_nama_q available for IT, by downloading everything, so
#     it is a bit slow (168 series)
# nm <- names(getSDMX('EUROSTAT', 'ei_nama_q.*.*.*.*.*.IT') )
# length(nm)  # 168

#  There are only Quarterly series in above, so next is the same (and also slow).
#  It works but several series have only  NaN values

#    tts <-  getSDMX('EUROSTAT', 'ei_nama_q.Q.*.*.*.*.IT') 
#    names(tts)

# for (i in 1: length(tts)) print( any(! is.nan(tts[[i]])))
# for (i in 1: length(tts)) print( sum(! is.nan(tts[[i]])))

# "ei_nama_q.Q.MIO-EUR.NSA.CLV2000.NA-B11.IT" %in% nm

# Retrieves but values are NaN
# tts2 <- getSDMX('EUROSTAT', "ei_nama_q.Q.MIO-EUR.NSA.CLV2000.NA-B11.IT") 

# This works and the series has data starting 1990Q1 (NaN prior to 1990)
#  tts2 <- getSDMX('EUROSTAT', "ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT") 


#  sdmxHelp()

getProviders()
#[1] "BIS"      "ILO"      "ECB"      "OECD"     "EUROSTAT"       with v1.1
#[1] "ILO"      "ECB"      "INEGI"    "OECD"     "EUROSTAT" "IMF" with v1.2

#with v1.3
#[1] "ABS"  "WB"  "ILO" "ECB" "OECD_RESTR" "NBB" "INEGI" "OECD" 
#                                        "UIS"  "EUROSTAT"   "IMF" 

# github version Sept 15, 2015 Build ID: 20150915-1113:
# [1] "ABS"        "OECD_RESTR" "EUROSTAT"   "ISTAT"      "INSEE"     
# [6] "WB"         "ILO"        "ECB"        "NBB"        "OECD"      
#[11] "INEGI"      "UIS"        "IMF"  

  # z <- getSDMX('EUROSTAT', 'ei_nama_q.Q.MIO-EUR.NSA.CLV2000.*.IT')[[1]]

 
############################ "BIS" ############################
# need account, not available to the publicly as of Nov 2014
#addProvider(name='BIS', endpoint='xxx', TRUE)


###########################################################################
###########################################################################

#  Notes regarding not yet providers 
#  These have SDMX but not sure about REST 

#  see also organizations listed at
# http://sdmx.org/wp-content/uploads/2014/09/SWG_members_8-9-2014.pdf

###########################################################################


######################################
#  RJSDMX function addProvider

######################################

## The addProvider function works only on SDMX 2.1 fully compliant providers. 
# All other versions of SDMX are "not so standard", and it is impossible (at 
# others are a 'custom' client

  addProvider(name='test', 
    endpoint='http://sdw-wsrest.ecb.europa.eu/service', FALSE)

  getFlows('test')
  
  z <- getSDMX('test', 'EXR.A.USD.EUR.SP00.A')


############ Swiss Federal Statistical Office ###############
#Site to the following e-mail address: webmaster@bfs.admin.ch
#http://www.bfs.admin.ch/bfs/portal/en/index/dienstleistungen/premiere_visite/07.html

# federal Swiss institute of statistics to disseminate their data in SDMX 
# and via a RESTful API.
# addProvider(name='SIS', 
#     endpoint='', FALSE)


############################ IStat ############################

# http://sodi.istat.it/sodiWS/service1.asmx.


############################ JEDH ############################

#http://www.jedh.org/jedh_dbase.html


############################ WHO ############################
# emailed gho_info @ who.int  Dec 7,2014

# their XLM looks like SDMX but does not say it is
# see
# http://apps.who.int/gho/data/node.resources
# http://apps.who.int/gho/data/node.resources.api?lang=en
# http://apps.who.int/gho/data/node.resources.examples?lang=en

# dimensions
#http://apps.who.int/gho/athena/api/ 

# example
# http://apps.who.int/gho/athena/api/GHO/WHOSIS_000001 
# http://apps.who.int/gho/athena/api/GHO/WHOSIS_000001?filter=COUNTRY:BWA 

################### Federal Reserve Board #####################

#Consumer credit from all sources (I think)
#https://www.federalreserve.gov/datadownload/Output.aspx?rel=G19&series=79d3b610380314397facd01b59b37659&lastObs=&from=01/01/1943&to=12/31/2010&filetype=sdmx&label=include&layout=seriescolumn


############################ Statistics Canada ############################


############################  Bank of Canada  ############################


############################  Fisheries and Oceans  ############################

# emailed for contact Dec 7, 1014

# FAO Fisheries has currently this SDMX 2.1 REST API with SDMX 2.0 messages:
# http://www.fao.org/figis/sdmx/
# FAO will publish this year also:
# http://data.fao.org/sdmx/

############################ British ONS  ############################

# https://www.ons.gov.uk/ons/apiservice/web/apiservice/home
# http://stackoverflow.com/questions/tagged/ons-api



############################ INE (Spain)  ############################
# reference here:
#http://www.bfs.admin.ch/bfs/portal/en/index/news/veranstaltungen/blank/blank/pax/04.parsys.1169.downloadList.27646.DownloadFile.tmp/countryreportspain.pdf

#Request for time series using the metadata
#description of the concepts that intervene
#(Population: dpop, Annual frequency
#FREQ.A and NUTS2 AREA2.ES11):
#http://servicios.ine.es/wstempus/SDMX/en/compact/dpop?metadata=FREQ.A:AREA_2.ES11: 

#addProvider(name='INE',
#    endpoint='http://servicios.ine.es/wstempus/', FALSE)

#addProvider(name='INE',
#    endpoint='http://servicios.ine.es/wstempus/SDMX/', FALSE)


#  getFlows('INE')

#sdmxHelp()

#tts <- getSDMX('INE', 'dpop')
#names(tts)

############################ UN ############################
# contaced on web form Dec 7, 2014
#http://unstats.un.org/unsd/tradekb/Knowledgebase/Data-Extraction-Using-Comtrade-Web-Service

##R example using json at
##http://comtrade.un.org/data/Doc/api/ex/r


# no SDMX yet but coming
# http://comtrade.un.org/data/doc/api/
# http://comtrade.un.org/data/doc/api/#Future


# UN Comtrade data request takes the following form:
# http://comtrade.un.org/api/get?parameters
# API call: 
#http://comtrade.un.org/api/get?max=50000&type=C&freq=A&px=HS&ps=2013&r=826&p=0&rg=all&cc=AG2&fmt=json

#http://comtrade.un.org/api/get?max=50000&type=C&freq=A&px=HS&ps=2013&r=826&p=0&rg=all&cc=AG2&fmt=sdmx

