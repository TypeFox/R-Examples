# ESRI Development Center (EDC) Entwickerforum Time & Space
# Muenster, 17. - 18.02.2011
# 
# Author: Daniel Nuest
# sos4R Webseite (downloads, news): http://www.nordholmen.net/sos4r
################################################################################

# SOS ##########################################################################
# Was muss ich wirklich wissen?
# Begriffe: offering, procedure, observedProperty, feature of interest
# Ablauf:	requestparameter sammeln, encoding, request senden, response parsen/
#	decoding.

# Workshop #####################################################################
# - Skript oeffnen in RGui, dann mit Strg + R aktuelle Zeile oder ausgewaehlten
# Code ausfuehren.
# - Fuer R-Neulinge gibt es Tipps und Tricks gleich mit :-).
# - Funktionen/Features und Baustellen von sos4R immer wieder zwischendurch.


##### Installation #############################################################
# Dependencies installieren (werden nicht aufgeloest wenn von Datei installiert 
# wird):
#install.packages("sos4R")
# Package laden
library("sos4R")
# sessionInfo()

# Weitere packages die benoetigt werden fuer Kartendarstellungen:
#install.packages("maps"); install.packages("mapdata"); install.packages("maptools")

##### Unterstuetzte und Standard-Features #######################################
SosSupportedOperations() # jeweils eine entsprechende R Funktion
SosSupportedConnectionMethods()

SosSupportedResponseFormats()
SosSupportedServiceVersions()
# Exchangeability feature: Parsing und encoding-Functionen koennen einfach 
# ausgetauscht werden (Liste mit Funktionen)

SosSupportedResponseModes()

SosSupportedResultModels()
## Exkurs: "O&M Profil": sos4R unterstuetzt, aus offensichtlichen Gruenden
# (52North Innovation Price for Geoinformatics), aber auch aus Ermangelung von
# allgemeinen Alternativen, das Profil des 52N SOS. Dies gilt auch fuer die 
# folgenden Operatoren:
SosSupportedSpatialOperators()
SosSupportedTemporalOperators()

SosDefaultConnectionMethod()
# SosDataFieldConvertingFunctions()
names(SosDataFieldConvertingFunctions())
SosDefaults()


##### Verbindung zu einem SOS erstellen ########################################
# Minimal-Server fuer Niederlande, nur 1 Monat Daten:
aqe <- SOS(url = "http://v-sos.uni-muenster.de:8080/SosAirQuality/sos")

# TIPP: Methoden beginnen mit 'sos...'
# > sos "TAB TAB" in Konsole
# sos "CTRL Space" in StatET
sosUrl(aqe)
sosVersion(aqe)
sosMethod(aqe)

# Standard R Funktionen werden oft unterstuetzt:
print(aqe)
summary(aqe)


##### Capabilities #############################################################
# Oeffne Seite im Browser: http://v-sos.uni-muenster.de:8080/SosAirQuality/sos?service=SOS&request=GetCapabilities
# Capabilities wurden jedoch bereits runtergeladen und geparsed:
sosCaps(aqe)

# OWS Common = OGC Web Services Common wurden implementiert fuer GetCapabilities
# Request und Exception Handling (OwsExceptionReport), modelliert als Klassen
# und parsing-Funktionen.

# TIPP: str() Funktion, siehe ?str, insbesondere mit max.level
str(sosCaps(aqe), max.level = 2)

# Teile der Capabilities:
sosContents(aqe)

#************#
# Aufgabe 01 #
#************#
# Wie koennen die weiteren Sektionen der Capabilities abgefragt werden?
# Die Elemente sind ServiceIdentification, ServiceProvider, OperationsMetadata,
# Filter_Capabilities.

sosFilter_Capabilities(aqe)
sosServiceIdentification(aqe)
sosServiceProvider(aqe)
sosOperationsMetadata(aqe)
str(sosOperationsMetadata(aqe))
sosOperation(aqe, "GetObservation"); sosOperations(aqe)[["GetObservation"]]

# Welche Information sind dort zu finden?
# Wie viel kostet die Benutzung diese SOS? --> sosServiceIdentification(...)
# Wer ist verantwortliche Kontaktperson? --> sosServiceProvider(...)


##### Zugriffsfunktionen #######################################################
# Weitere wichtige Elemente der Capabilities, sind ueber getter bzw.  accessor 
# abrufbar, da sie spaeter zum erstellen der Requests benoetigt werden:
sosOfferingIds(aqe) # == names(sosOfferings(aqe))

sosOfferings(aqe) # > print(sosOfferings(aqe))

# Indexierungsarten (generisch R):
sosOfferings(aqe)[2:3]				# Subliste
sosOfferings(aqe)[[3]]				# einzelnes Element
sosOfferings(aqe)[["PM10"]]		# ueber die ID 

pm10.off <- sosOfferings(aqe)[["PM10"]]

sosId(pm10.off)
sosProcedures(pm10.off)
sosObservedProperties(pm10.off)
sosBoundedBy(pm10.off, bbox = TRUE)
sosTime(pm10.off)
sosGetCRS(pm10.off)


#************#
# Aufgabe 02 #
#************#
# Welche der Funktionen aus diesem Abschnitt funktionieren auch fuer SOS-Objekte?
# Was gibt es dabei zu beachten?

sosTime(aqe) # benutzt GetObservation operation description
sosObservedProperties(aqe)
str(sosObservedProperties(aqe))
sosProcedures(aqe)
sosGetCRS(aqe)

#sosId(aqe) # GEHT NICHT! SOS hat keine ID.
#sosBoundedBy(aqe) # GEHT NICHT! BoundingBox ist Offering-spezifisch.

# Welches der Offerings hat die meisten Sensoren?
# TIPP: ?length

# Loesung 1:
lapply(X = sosProcedures(aqe), FUN = length)
lapply(X = sosProcedures(sosOfferings(aqe)), FUN = length)

# Loesung 2:
procedures <- sosProcedures(sosOfferings(aqe))
length(procedures[[1]])
length(procedures[[2]])
length(procedures[[3]])

# Loesung 3:
length(sosProcedures(sosOfferings(aqe)[[1]]))
length(sosProcedures(sosOfferings(aqe)[[2]]))
length(sosProcedures(sosOfferings(aqe)[[3]]))


##### Grafische Ausgabe von SOS und Offerings ##################################
plot(pm10.off) # alleinstehend nicht sehr hilfreich

# Hintergrundkarte erstellen -- muss nicht verstanden werden:
library(maps); library(mapdata); library(maptools)
data(worldHiresMapEnv)
crs <- unique(sosGetCRS(aqe))[[1]]
region <- map.where(database = "world",
		sosCoordinates(pm10.off)) # Finde region basierend auf einem offering.
map <- pruneMap(map(database = "world", region = region,
				plot = FALSE))
map.lines <- map2SpatialLines(map, proj4string = crs)

plot(map.lines, col = "grey50", lwd = 1)
plot(pm10.off, add = TRUE, lwd = 3)
title(main = paste("Offering '", sosId(pm10.off), "' at", sosTitle(aqe),
				sep = ""),
		sub = paste(sosAbstract(aqe), " --- ",
				toString(names(sosOfferings(aqe)))))

#************#
# Aufgabe 03 #
#************#
# Wie koennen alle Offerings geplottet werden?
# TIPP fuer Fortgeschrittene: ?text

# Loesung 1: "copy und paste"
plot(map.lines, col = "grey50", lwd = 1)
plot(sosOfferings(aqe)[[1]], add = TRUE, lwd = 3, border = "green")
plot(sosOfferings(aqe)[[2]], add = TRUE, lwd = 3, border = "blue")
plot(sosOfferings(aqe)[[3]], add = TRUE, lwd = 3, border = "red")

# Loesung 2:
plot(map.lines, col = "grey50", lwd = 1) # neuen Hintergrundplot erzeugen.
plot(aqe, add = TRUE, lwd = 2)


##### Sensor Metadaten abfragen ################################################
sensor2 <- describeSensor(aqe, sosProcedures(pm10.off)[[2]])
sensor2		
sensor2@xml # Viel Information, einzelne Details ueber Getter abfragbar, dies
			# verlangt aber Konformitaet mit SensorML Profile for Discovery

# Getter:
sosId(sensor2)
sosBoundedBy(sensor2)

#************#
# Aufgabe 04 #
#************#
# Wie ist die Bounding Box von 'sensor2'?

sosBoundedBy(sensor2)

# Wo in Deutschland (Koordinaten und/oder Plot) ist 'sensor2'?

sosCoordinates(sensor2)
# geht davon aus dass vorheriger plot nicht geschlossen wurde
plot(sensor2, add = TRUE, pch = 7, cex = 2)


##### Messungen abfragen #######################################################
# TIPP: "inspect" benutzen um die requests und responses kennen zu lernen!

##### Temporaere Ausschnitte ####################################################
# 1 Stunde im Dezember, alle Daten zu viel wegen parallelen Anfragen!
# TIPP: sosCreate... - Funktionen helfen korrekte Strukturen zu erzeugen:
aug2007.6Hrs = sosCreateEventTimeList(sosCreateTimePeriod(sos = aqe,
				begin = as.POSIXct("2007/08/10 08:00"),
				end = as.POSIXct("2007/08/10 14:00")))

# Daten abrufen, nur zeitlicher Auschnitt:
aug2007.obs <- getObservation(sos = aqe, inspect = TRUE,
		offering = pm10.off, eventTime = aug2007.6Hrs)
warnings()

# PROBLEM: Nicht unterstuetzte Datenfelder werden abgefragt! Konverter ergaenzen:
aqe.converters <- SosDataFieldConvertingFunctions(
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[PM10]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration.NO2." = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[O3]" = sosConvertDouble,
		"http://www.opengis.net/def/property/OGC/0/SamplingTime" = sosConvertTime,
		"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString)

# Neue Verbindung erstellen, dieses mal mit Konvertern:
aqe <- SOS(sosUrl(aqe), dataFieldConverters = aqe.converters)

# Nochmaliges daten abrufen, diesmal alternativ auf Basis der offering ID:
# Daten fuer ganz Deutschland!
aug2007.obs <- getObservation(sos = aqe, offering = "NO2", #verbose = TRUE,
		eventTime = aug2007.6Hrs)
# Parsing dauert laenger als request -> wenig Daten in vielen Observations, hoher
# XML overhead?

# Nur einige Stationen, mit 'inspect' Parameter
aug2007.obs.small <- getObservation(sos = aqe, offering = "NO2",# verbose = TRUE,
		inspect = TRUE,
		procedure = sosProcedures(aqe)[[1]][10:20],
		eventTime = aug2007.6Hrs)

##### Response erforschen ######################################################
print(aug2007.obs.small)
str(aug2007.obs.small, max.level = 5) # limitieren der Strukturtiefe

## OmObservationCollection und OmObservation modellieren O&M Datenmodell
# generisch in R. Dies gibt es in begrenztem Umfang auch fuer GML, SWE, SensorML,
# OGC, SA, ... Namespaces.
## Mittelfristig denkbar: Nicht nur De- und Encoding der Elemente wie sie fuer 
# den SOS benoetigt werden, sondern generische Transformationsfunktionen
# (Stichwort "coercion") zwischen O&M/GML/SWE Datenmodell und sp, base-R, ...

# ObservationCollection ist auf vielfaeltige Art indizierbar:
aug2007.obs[[1]]
aug2007.obs[10:14]
aug2007.obs[c(17,21)]

# Koordinaten, Features und BoundingBox abfragen:
sosCoordinates(aug2007.obs[[1]])
sosCoordinates(aug2007.obs[10:20])
sosFeatureIds(aug2007.obs)[10:14]
sosBoundedBy(aug2007.obs, bbox = TRUE)

#************#
# Aufgabe 05 #
#************#
# Welche procedures, observedProperties und features sind in den erhaltenen
# Observations zu finden? TIPP: Zugriffsfunktionen!

sosProcedures(aug2007.obs)
sosObservedProperties(aug2007.obs[[10]])

# Wie koennen die Koordinaten fuer die Observations 10 bis 20 abgefragt werden?

sosCoordinates(aug2007.obs[10:20])

# Was ist der Unterschied zwischen sosFeatureIds(aug2007.obs)[10:12] und 
# sosFeatureIds(aug2007.obs[10:12])?

sosFeatureIds(aug2007.obs)[10:12]
# Fragt alle feature ids ab und nimmt von dieser zusammengefuegten Liste die
# Elemente 42 bis 44:
str(sosFeatureIds(aug2007.obs)[10:12]) # List of 3 character

sosFeatureIds(aug2007.obs[10:12])
# Fragt observations 10 bis 14 ab und ruft auf dieser Liste (!) fuer jedes
# Element die Funktion sosFeatureIds(...) auf.
str(sosFeatureIds(aug2007.obs[10:12])) # List of 3 Lists


##### Tatsaechliche Daten erforschen ############################################
# sosResult(...) ist die wichtigste Methode:
sosResult(aug2007.obs[[1]])

# Daten enthalten Zeitreihe an jedem Punkt.

# Testplot
plot(sosResult(aug2007.obs[[1]]))

# Was ist die Struktur des Results?
class(sosResult(aug2007.obs[[1]]))

# R Hilfe zu data.frame (und allen anderen Funktionen)
?data.frame
help("data.frame")
# HTML help:
help.start()

# Wie kann ich die Daten verschiedene Stationen zusammenfuegen?
sosResult(aug2007.obs[20:21])			# Result der 20. und 21. Observations
sosResult(aug2007.obs)[20:21,]		# Zeile 20, 21 von allen (!) Daten

aug2007.result <- sosResult(aug2007.obs)

# Warum XML parsen, wenn sowieso nur CSV-Werte zur einer "Tabelle" geparst
# werden? Attribute enthalten Metadaten aus der Observation!
names(aug2007.result)
attributes(aug2007.result[["Concentration.NO2."]])

# Kombination der results mit den Koordinaten
aug2007.coords <- sosCoordinates(aug2007.obs)[1:5,]
merge(x = aug2007.result, y = aug2007.coords)[1:3,]

# Geht viel einfacher!
aug2007.data <- sosResult(aug2007.obs, coordinates = TRUE)
aug2007.data[1:3,]

#************#
# Aufgabe 06 #
#************#
# Was sind der maximale/minimale, der Durchschnittswert, der Median und die
# Quantile von NO2 fuer alle heruntergeladenen Daten?

# Loesung 1:
max(aug2007.data[["Concentration.NO2."]])
min(aug2007.data[["Concentration.NO2."]])
mean(aug2007.data[["Concentration.NO2."]])
median(aug2007.data[["Concentration.NO2."]])
quantile(aug2007.data[["Concentration.NO2."]])

# Loesung 2:
summary(aug2007.data)
summary(aug2007.data[["Concentration.NO2."]])

# Wie ist das Phaenomen Concentration.NO2. definiert?

# Antwort: "Amount of nitrogen dioxide (NO2) as a fraction of host medium"

# Loesung 1:
# Browser: http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration.NO2.

# Loesung 2: XML parsing
# Nur um zu zeigen dass Erweitern mit Paket XML nicht so schlimm ist:
definition <- getURL("http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration.NO2.")
definition.xml <- xmlParse(definition)
getNodeSet(doc = definition.xml, path = "//gml:description/text()")[[1]]

##### Thematische Ausschnitte ##################################################

# Ein Phaenomen (auch wenn es sowieso nur eines ist) fuer eine Station, ein Jahr:
no2.off <- sosOfferings(aqe)[["NO2"]]
time.2007 = sosCreateEventTimeList(sosCreateTimePeriod(sos = aqe,
				begin = as.POSIXct("2007/01/01"),
				end = as.POSIXct("2007/12/31")))

sosProcedures(no2.off)
# Beachte das Benennungsschema: "DE" "NW", "NL"
# Stationsuebersicht: http://www.eea.europa.eu/themes/air/airbase/map-stations
# Station auswaehlen -> Meta Data
myStationID <- "NL00639" # "DESH024"
idx <- grep(pattern = myStationID, x = sosProcedures(no2.off))
myStation <- sosProcedures(no2.off)[idx]
myStation

sosObservedProperties(no2.off)

# Filtern mit observed property nicht notwendig, da sowieso ein Offering nur 
# ein Phaenomen liefert.
obs.myStation.2007 <- getObservation(sos = aqe,
		offering = no2.off, # inspect = TRUE,
		procedure = myStation,
		observedProperty = sosObservedProperties(no2.off),
		eventTime = time.2007)

result.myStation.2007 <- sosResult(obs.myStation.2007)
summary(result.myStation.2007)

# Zeitreihenplot:
plot(result.myStation.2007[["SamplingTime"]],
		result.myStation.2007[["Concentration.NO2."]],
		type = "l",
		main = paste("NO2 in", myStation,
				min(result.myStation.2007[["SamplingTime"]]), "-",
				max(result.myStation.2007[["SamplingTime"]])),
		sub = myStation,
		xlab = "Time",
		ylab = paste("NO2 (",
				attributes(result.myStation.2007[["Concentration.NO2."]])[["unit of measurement"]],
				")", sep = ""))
locRegr <- loess(result.myStation.2007[["Concentration.NO2."]]~
				as.numeric(result.myStation.2007[["SamplingTime"]]),
		result.myStation.2007, enp.target = 5)
p = predict(locRegr)
lines(p ~ result.myStation.2007[["SamplingTime"]], col = 'blue',lwd = 4)

# Histogramm:
hist(result.myStation.2007[["Concentration.NO2."]])


##### Raeumliche Ausschnitte ####################################################
#SosSupportedSpatialOperators()
sosBoundedBy(no2.off)

# NRW
#request.bbox <- sosCreateBBOX(lowLat = 49.84, lowLon = 5.98,
#		uppLat = 52.12, uppLon = 10.15, srsName = "urn:ogc:def:crs:EPSG:4326")
# Region Amsterdam:
request.bbox <- sosCreateBBOX(lowLat = 52.276, lowLon = 4.667,
		uppLat = 52.450, uppLon = 5.049, srsName = "urn:ogc:def:crs:EPSG::4326")
request.bbox.foi <- sosCreateFeatureOfInterest(spatialOps = request.bbox)
request.bbox.foi

# Alle NO2 Daten von 2007 in Bounding Box
obs.2007.bbox <- getObservation(sos = aqe, # inspect = TRUE,
		offering = no2.off,
		featureOfInterest = request.bbox.foi,
		eventTime = sosCreateEventTimeList(sosCreateTimePeriod(sos = aqe,
						begin = as.POSIXct("2007/01/01"),
						end = as.POSIXct("2007/12/31"))))
obs.2007.bbox
sosBoundedBy(obs.2007.bbox, bbox = TRUE)

#************#
# Aufgabe 07 #
#************#
# Wann und wo (Koordinaten) sind Daten des Offerings NO2 im Vergleich zu den
# abgefragten Daten verfuegbar?

# Wo:
sosBoundedBy(no2.off, bbox = TRUE)
summary(sosCoordinates(obs.2007.bbox)[c("lat","lon")])

# Wann:
result.bbox <- sosResult(obs.2007.bbox)
range(result.bbox[["SamplingTime"]])
sosTime(no2.off)

# Wie viele Messstationen gibt es in der bounding box, die NO2-Werte liefern?
length(sosProcedures(obs.2007.bbox))

##### Daten -> sp ##############################################################
# sp-Object bzw. -Klassen sind die komfortablen Schnittstellen (im Vergleich zu
# einfachen data.frames) von sos4R in andere Pakete zur Raum-zeitlichen Analyse.

result.bbox <- sosResult(obs.2007.bbox, coordinates = TRUE)
obs.crs <- sosGetCRS(obs.2007.bbox)

# Die Spalten lon, lat werden fuer die Koordinaten des SPDF, die anderen Spalten
# fuer die Daten des SPDF benutzt.
no2.spdf.bbox <- SpatialPointsDataFrame(
		coords = result.bbox[,c("lon", "lat")],
		data = result.bbox[,
				c("SamplingTime", "feature", "Concentration.NO2.")],
		proj4string = obs.crs)
summary(no2.spdf.bbox)

# Viele Funktionen aus sp, ... nun verfuegbar
bbox(no2.spdf.bbox)

# Im Vergleich zu den Daten ohne bounding box, hierbei wird eine "Abkuerzung"
# verwendet wenn die Spaltennamen bekannt sind (sind sie bei 52N SOS)
obs.no2.2007 <- getObservation(sos = aqe, # inspect = TRUE,
		offering = no2.off)
result.no2.2007 <- sosResult(obs.no2.2007, coordinates = TRUE)
no2_spdf <- as(obs.no2.2007, "SpatialPointsDataFrame")
summary(no2_spdf)
bbox(no2_spdf)

# Coercion einer einzelnen Observation ist ebenfalls moeglich, jedoch bis jetzt
# auf SpatialPointsDataFrame beschraenkt. Hier ist natuerlich aufwendiges/
# komfortables/maechtiges Mapping zwischen O&M Datenstrukturen und sp/R/spacetime
# Datenmodellen ein langfristiges Ziel. Potentiell sogar in beide Richtungen,
# (Transactional SOS Profile). Funktioniert (nur) fuer 52N Profil oder aehnliche.
spdf.1 <- as(obs.2007.bbox[[1]], "SpatialPointsDataFrame")
summary(spdf.1)
levels(spdf.1[["FeatureOfInterest"]])


#************#
# Aufgabe 08 #
#************#
# Wo sind die Messtationen von no2.spdf?

coordinates(no2_spdf) # nur Koordinates, coordinates ist eine sp-Funktion
plot(x = map.lines, col = "grey")
plot(no2_spdf, pch = 20, col = "blue", add = TRUE)

# Frage Daten fuer eine beliebige Woche ab und erzeuge einen data.frame, benutze
# auch den SOS fuer Deutschland. Wichtig: kleiner zeitlicher oder raeumlicher
# Ausschnitt, damit der Service nicht ueberansprucht wird.
aqe.de <- SOS(url = "http://giv-uw.uni-muenster.de:8080/AQE/sos",
		dataFieldConverters = SosDataFieldConvertingFunctions(aqe.converters))

# ... viel Spass beim Programmieren!

##### Demos ####################################################################
demo(package = "sos4R")

# Demos laufen lassen (enhalten weiterfuehrende Beispiele mit plots usw.):
#demo("weathersos")
#demo("pegel")
#demo("airquality")


##### Fragen? ##################################################################
vignette("sos4R")
sosCheatSheet()
# Mailingliste: http://52north.org/resources/mailing-list-and-forums/
# Forum:		http://geostatistics.forum.52north.org/
# Webseite:		http://www.nordholmen.net/sos4r
# Kontakt:		d.nuest@52north.org
