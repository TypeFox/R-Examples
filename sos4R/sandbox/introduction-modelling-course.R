# V+ue Einfuehrung Modellierung dynamischer raeumlicher Prozesse, WS 2010/2011
# 
# Author: Daniel Nuest
###############################################################################

# Was ist OGC?
#	- WMS, KML, GML, WFS, WCS, CSW, OwsCommon
# Was ist ein SOS?
# - OGC Standard fuer Onlinespeicherung von Sensormessungen:
#	http://www.opengeospatial.org/standards/sos
#	- Gute uebersicht: http://de.wikipedia.org/wiki/Sensor_Observation_Service
#	- Begriffe: FOI, Observation, Offering, Phenomenon, Procedure, In-Situ
#	- Verbindungsarten: HTTP GET, HTTP POST, SOAP
#	- Dokumentenkodierung: SensorML, O&M, SweCommon, ...

# Warum muss ich mich damit beschaeftigen?
# - SWSL: http://swsl.uni-muenster.de
# - Standards = Interoperabilitaet
#	(semantische ~ --> http://musil.uni-muenster.de/)

# Was muss ich wirklich wissen?
# Begriffe:
#	- offering
#	- procedure
#	- observedProperty
#	- feature of interest
# - Aber im Endeffekt ist es doch von Fall zu Fall unterschiedlich.

# Kann ich mit R nicht Daten einfach so runterladen?
library("RCurl")
temp <- getURL("ifgi.uni-muenster.de")
cat(temp)
temp <- getURL("http://giv-sos.uni-muenster.de:8080/52nSOSv3/sos?request=GetCapabilities&version=1.0.0&service=SOS")
print(temp)
cat(temp)
# CSV Datei runterladen und importieren
download.file(url = "http://geography.uoregon.edu/geogr/data/csv/cities.csv",
		destfile = "cities.csv")
cities <- read.csv(paste(getwd(), "/cities.csv", sep = ""))
# Eine einfache Tabelle - kein Problem!

#  10080 000     SVALBARD LUHTHAVN                   78.3n   15.5e    29m  temperature in degrees centigrade (  99.0 is missing)   
#  10080 001  par    jan    feb    mar    apr    may    jun    jul    aug    sep    oct    nov    dec    djf    mam    jja    son    ann    
#  10080 1977   3   99.0   99.0   99.0   99.0   99.0   99.0   99.0   99.0   99.0   -4.1   -6.9   -9.5   99.0   99.0   99.0   99.0   99.0
#  10080 1978   3  -20.7  -20.4  -17.7   -9.9   -3.7    1.7    6.0    5.0   -1.8   -5.8  -10.9  -12.4  -16.9  -10.4    4.2   -6.2   -7.6
#  10080 1979   3  -20.5  -20.6  -16.0  -17.4   -8.0   -0.5   99.0    4.8    1.0   -4.3   99.0  -12.7  -17.8  -13.8   99.0   99.0   99.0
# [...]  10080 2003   3  -18.6  -12.3  -17.3   -9.6   -2.5    2.8    7.0    6.4    0.6   -4.3   -6.7  -19.1  -12.8   -9.8    5.4   -3.5   -6.1
#  10080 2003   3  -18.6  -12.3  -17.3   -9.6   -2.5    2.8    7.0    6.4    0.6   -4.3   -6.7  -19.1  -12.8   -9.8    5.4   -3.5   -6.1
#  10080 2004   3  -17.0  -17.5   -7.7   -4.3   -1.8    3.0    7.6    5.7    1.3   -3.6  -11.5   -6.6  -17.9   -4.6    5.4   -4.6   -4.4
#  10080 991 #yrs    27.    27.    27.    27.    27.    27.    26.    27.    27.    28.    27.    28.    27.    27.    26.    26.    26.
#  10080 992 mean  -14.7  -15.0  -13.6  -10.7   -3.5    2.5    6.3    5.2    0.8   -4.7   -9.1  -12.5  -14.2   -9.3    4.7   -4.4   -5.7
#  10080 993 stdv    4.4    3.8    3.7    3.3    1.3    1.0    0.9    0.9    1.7    3.1    3.8    4.2    3.0    2.1    0.7    2.0    1.3
# Format?!?
# http://dss.ucar.edu/datasets/common/wmssc/format.html

# --> Vorteile vom SOS ausnutzen: Wohldefiniertes Markup, Queries

# Was ist sos4R?
# - Webseite (downloads, news): http://www.nordholmen.net/sos4r
# - Funktionalitaet
# 	- Core Profile (plug GetObservationById)
#		- GetCapabilities
#		- GetObservation
#		- DescribeSensor
#	- GET und POST
#	- Austauschbarkeit

##### Installation #############################################################
# Dependencies installieren (werden nicht aufgeloest wenn von Datei installiert 
# wird):
install.packages("RCurl")
install.packages("XML")

# Package herunterladen
pkgName = "sos4R_0.1-08.tar.gz"
download.file(url = paste(
				"http://www.nordholmen.net/sos4r/download/", pkgName, sep =""),
		destfile = pkgName)
# Package installieren
install.packages(paste(getwd(), "/", pkgName, sep = ""))

# Package laden
library("sos4R")
# sessionInfo()

# Unterstuetzte Features (werden hoffentlich in Zukunft mehr)
SosSupportedOperations()
SosSupportedServiceVersions()
SosSupportedConnectionMethods()
SosSupportedResponseFormats()
SosSupportedResponseModes()
SosSupportedResultModels()
SosSupportedSpatialOperators()
SosSupportedTemporalOperators()

# Default Features
SosDefaultConnectionMethod()
SosDataFieldConvertingFunctions()
SosDefaults() # TODO named list with defaults


##### Verbindung zu einem SOS erstellen ########################################
mySOS = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")

# TIPP: Methoden beginnen mit 'sos...'
# > sos "TAB TAB" in Konsole
# sos "CTRL Space" in StatET
sosUrl(mySOS)
sosVersion(mySOS)
sosTimeFormat(mySOS)
sosMethod(mySOS)

##### Capabilities abfragen ####################################################
# http://v-swe.uni-muenster.de:8080/WeatherSOS/sos?service=SOS&request=GetCapabilities
# Wurden bereits runtergeladen:
sosCaps(mySOS)
# Originaldokument:
sosCapabilitiesDocumentOriginal(mySOS)

# Capabilities Erforschen:
sosContents(mySOS)
sosFilter_Capabilities(mySOS)
sosServiceIdentification(mySOS)
sosServiceProvider(mySOS) # @serviceContact

# Wichtig:
sosOfferings(mySOS)
off.temp <- sosOfferings(mySOS)[["ATMOSPHERIC_TEMPERATURE"]]
sosOfferingIds(mySOS)
names(sosOfferings(mySOS))

sosId(off.temp)
sosOfferings(mySOS)[1:3]

sosProcedures(mySOS)
sosProcedures(off.temp)

sosObservedProperties(mySOS)
sosObservedProperties(off.temp)

sosBoundedBy(off.temp)
str(sosBoundedBy(off.temp)) # Nicht so schoen ...

sosTime(mySOS)
off.temp.time <- sosTime(off.temp)
str(off.temp.time) # modelliert XML
# "wirklichen" Startzeitpunkt abfragen
off.temp.time@beginPosition@time
class(off.temp.time@beginPosition@time)

##### Sensor Metadaten abfragen ################################################
describeSensor(sos = mySOS, procedure = "meinTollerSensor")
# Warnung! Und: Antwort ist OwsExceptionReport!

# TIPP: verbose option in service operation functions
describeSensor(sos = mySOS, procedure = "lala", verbose = TRUE)

sensor2 <- describeSensor(mySOS, sosProcedures(off.temp)[[2]])
sensor2
str(sensor2)
sensor2@xml

##### Messungen abfragen #######################################################
# Einfachster Fall: "latest observation" fuer ein gesamtes offerings
obs.temp.latest <- getObservation(sos = mySOS, offering = off.temp,
		latest = TRUE, inspect = TRUE)
# TIPP: "inspect" benutzen um die requests und responses kennen zu lernen!

##### Response erforschen ######################################################
# print Methoden
obs.temp.latest
# TIPP: str(...) fuer Einblick unter die Motorhaube

# ObservationCollection behaves like a list in most cases
length(obs.temp.latest)
obs.temp.latest[[1]]
obs.temp.latest[2:5]

# Koordinaten, Features und BoundingBox abfragen
sosCoordinates(obs.temp.latest)
sosCoordinates(obs.temp.latest[[1]])
sosFeatureIds(obs.temp.latest)
sosBoundedBy(obs.temp.latest)


##### Daten erforschen #########################################################
# sosResult(...) ist die wichtigste Methode
sosResult(obs.temp.latest[[2]])
obs.temp.latest.result <- sosResult(obs.temp.latest[1:2])

# Nur ein ganz normaler data.frame ... Attribute enthalten Metadaten. Diese 
# gehen nach dem "merge" verloren!
attributes(obs.temp.latest.result[["urn:ogc:def:property:OGC::Temperature"]])

# Kombination der results mit den Koordinaten
obs.temp.latest.coords <- sosCoordinates(obs.temp.latest)
obs.temp.latest.data <- merge(x = obs.temp.latest.result,
		y = obs.temp.latest.coords)
obs.temp.latest.data

##### Temporaere Ausschnitte ####################################################
# Erstellen der event time fuer GetObservation requests:
period.august09 <- sosCreateEventTimeList(
		sosCreateTimePeriod(sos = mySOS,
				begin = as.POSIXct("2009-08-01 00:00"),
				end = as.POSIXct("2009-08-31 12:00")))
#sosCreateTimeInstant()
#sosCreateTimePeriod()
str(period.august09)
#SosSupportedTemporalOperators()

obs.august09 <- getObservation(sos = mySOS,
		offering = off.temp,
		procedure = sosProcedures(off.temp),
		eventTime = period.august09)

obs.temp.august09.result <- sosResult(obs.august09)
summary(obs.temp.august09.result)
str(obs.temp.august09.result)
obs.temp.august09.result[100:103,]

##### Raeumliche Ausschnitte ####################################################
#SosSupportedSpatialOperators()

sosBoundedBy(off.temp)
#request.bbox <- sosCreateBBOX(lowLat = 10.0, lowLon = 2.0,
#		uppLat = 30.0, uppLon = 5.0, srsName = "urn:ogc:def:crs:EPSG:4326")
request.bbox <- sosCreateBBOX(lowLat = 50.0, lowLon = 7.0,
		uppLat = 52.0, uppLon = 9.0, srsName = "urn:ogc:def:crs:EPSG:4326")
request.bbox.foi <- sosCreateFeatureOfInterest(spatialOps = request.bbox)
request.bbox.foi

obs.august09.bbox <- getObservation(sos = mySOS,
		offering = off.temp,
		featureOfInterest = request.bbox.foi,
		eventTime = period.august09,
		inspect = TRUE)
obs.august09.bbox
sosCoordinates(obs.august09.bbox)

# Beliebige raeumliche Filter ueber FOI sind auch moeglich
# -> siehe SOS Spec.
# -> siehe sos4R Code

##### Fortgeschrittenes Filtern ################################################
# Mit feature of interest:
off.temp.fois <- sosFeaturesOfInterest(off.temp)
request.fois <- sosCreateFeatureOfInterest(objectIDs = list(off.temp.fois[[1]]))

obs.august09.fois <- getObservation(sos = mySOS, offering = off.temp,
		featureOfInterest = request.fois,
		eventTime = period.august09)
obs.august09.fois
# Object of class OmObservationCollection with 1 members. 

# Weitere Filter sind derzeit nicht implementiert...

##### Neue Konverter ###########################################################
# GET Verbindung
MBARI <- SOS("http://mmisw.org/oostethys/sos",
		method = SosSupportedConnectionMethods()[["GET"]])
myOff <- sosOfferings(MBARI)[[1]]
myProc <- sosProcedures(MBARI)[[1]]
mbariObs <- getObservation(sos = MBARI, offering = myOff, procedure = myProc,
		inspect = TRUE)
# Warnmeldungen!

?SosDataFieldConvertingFunctions

# So geht es:
myConverters <- SosDataFieldConvertingFunctions(
		# mapping for UOM:
		"C" = sosConvertDouble,
		"S/m" = sosConvertDouble,
		# mapping for definition:
		"http://mmisw.org/ont/cf/parameter/sea_water_salinity" = sosConvertDouble)
MBARI <- SOS("http://mmisw.org/oostethys/sos",
		method = SosSupportedConnectionMethods()[["GET"]],
		dataFieldConverters = myConverters)
myOff <- sosOfferings(MBARI)[[1]]
myProc <- sosProcedures(MBARI)[[1]]
mbariObs <- getObservation(sos = MBARI, offering = myOff, procedure = myProc)

sosResult(mbariObs)

##### Daten -> sp/gstat ########################################################

# Typisch fuer das Geoinformatikerleben: Anwenden von Software anderer Leute...
# Viel Erfolg!

##### Demos ####################################################################
demo(package = "sos4R")
# Demos laufen lassen (enhalten weiterfuehrende Beispiele mit plots usw.):
demo("weathersos")
demo("pegel")		# Moeglicherweise SOS defekt/offline/langsam
demo("airquality")	# Moeglicherweise SOS defekt/offline/langsam
demo("marinemeta")


##### Fragen? ##################################################################
# Mailingliste: http://52north.org/resources/mailing-list-and-forums/
# 				Bitte den Posting-Guide beachten!
# Forum: http://geostatistics.forum.52north.org/
# Kontakt: daniel.nuest@uni-muenster.de

##### Mach mit! ################################################################
# Entwicklerdokumentation und To-Do-Liste:
# https://wiki.52north.org/bin/view/Geostatistics/Sos4R#To_Do_List
