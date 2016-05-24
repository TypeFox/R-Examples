exampleCommands <- function(){
	doitall <- TRUE
	#doitall <- FALSE
	if(doitall){

###ADD DATA VALUES
	try(getMetadata("Site"), silent=TRUE)
	addSite(Code="testSpatialPoints", Name="Virtual test site", x=25, y=56,
	LatLongDatum="WGS84", Elevation=350, State="Germany")
	addVariable(Name="Temperature, transducer signal", Unit="degree celsius", ValueType="Field Observation",
	GeneralCategory="Hydrology", Code="test_temp")
	addQualityControlLevel(ID=2,Code="test_ok", Definition="The default values")
	addISOMetadata(TopicCategory="Unknown", Title="Testdata",
		Abstract="This data is created to test the functions of RObsDat")
	addSource(Organization="Your Org", SourceDescription="Madeup data", 
		SourceLink="RObsDat Documentation", ContactName="Yourself",
		Metadata="Testdata")
	
	
	example.data <- xts(1:40, seq(as.POSIXct("2014-01-01", tz="UTC"), 
		as.POSIXct("2014-02-09", tz="UTC"), length.out=40))
	example.data[40] <- 30
	example.data[35] <- 22

	addDataValues(example.data[1:20], Site="testSpatialPoints", Variable="test_temp",  
	Source="Madeup", QualityControlLevel="test_ok")
	#Avoid duplicates autmatically
	addDataValues(example.data, Site="testSpatialPoints", Variable="test_temp",  
	Source="Madeup", QualityControlLevel="test_ok")
	
	
### TEST STPLOT
	testSiteData <- getDataValues(Site="test")
	stplot(testSiteData, mode="ts")

	
### TEST SELECTION
	#test selection Main Data
	selectedData1 = testSiteData[, 10:20]
	selectedData2 = testSiteData[,,'Temperatur, transducer signal']

	#test selection for other Slots
	selectedData3 = testSiteData@ValueIDs[,10:20]
	selectedData4 = testSiteData@ValueIDs[,,'Temperatur, transducer signal']

	
### VERSION MANAGEMENT

### UPDATE
	to.correct <- which(testSiteData@data > 30)
	testSiteData@data[to.correct,] <- 20
	testSiteData@data[39,] <- 32

	updateDataValues(testSiteData, "Correction of wrong value")

	ver2 <- testSiteData
	ver2@data[10:13,] <- 60
	updateDataValues(ver2, "Changing more data")

	ver3 <- testSiteData
	ver3@data[30:32,] <- 33
	updateDataValues(ver3, "Ups, I used 60 instead of 33 by mistake")

	
### DELETE
	# via ValueID:
	deleteDataValues(testSiteData@ValueIDs[,36],  "And finally remove one value via ID")
	deleteDataValues(testSiteData@ValueIDs[,20:26],  "And finally remove several values via ID")
	
	# via direct access
	deleteDataValues(testSiteData[,30:35],  "And finally remove multi values via ID")
	
	to.delete <- testSiteData@data == 60
	if(any(to.delete)){
		deleteDataValues(testSiteData@ValueIDs[,to.delete],  "And finally remove several value")
	}
	
### GET OLD VERSION
	getDataVersions()

	#show no deleted values
	versionQuery <- getDataValues(Site=1, VersionID=1)
	stplot(versionQuery, mode = 'ts')
	
	versionQuery <- getDataValues(Site=1, VersionID=2)
	stplot(versionQuery, mode = 'ts')
	
	#show deleted values
	versionQuery <- getDataValues(Site=1, VersionID=2, show.deleted=TRUE)
	stplot(versionQuery, mode = 'ts')
	}

}

cleanupMySQL <- function(con){
		for(i in 1:6){
			try(dbGetQuery(con, "DROP TABLE if exists DataValues, DataValuesRepository;"), silent=TRUE)
			try(dbGetQuery(con, "DROP TABLE if exists `Categories`,  `DerivedFrom`, `GeneralCategoryCV`, `GroupDescriptions`, `Groups`, `ISOMetadata`, `LabMethods`, `Methods`, `ODMVersion`, `OffsetTypes`, `Qualifiers`, `QualityControlLevels`, `SampleMediumCV`, `Samples`, `SampleTypeCV`, `SeriesCatalog`, `Sites`, `Sources`, `SpatialReferences`, `SpeciationCV`, `TopicCategoryCV`, `Units`, `ValueTypeCV`, `VariableNameCV`, `Variables`, `Versions`, `VerticalDatumCV`, CensorCodeCV, DataTypeCV , DataValues, Synonyms;") , silent=TRUE)
		}
}



longExample <- function() {
options(error = recover)
	if(!as.logical(getOption("testLongExample", FALSE))){
		cat("Not testing long example. Use 'options(testLongExample = TRUE)' to do so\n")
		return()
	}
	# prepare countries --> only once!
	options(geonamesUsername="cite.all")
	if(!requireNamespace("geonames", quietly=TRUE)){
		cat("Please install geonames with install.packages(geonames) to run this test\n")
		return()
	}
	cat("Importing Geonames from web\n")
	all.countries <- geonames::GNcountryInfo(lang="EN")

	cat("Importing other Geonames\n")
	addSite(Code=all.countries$isoAlpha3, Name=all.countries$countryName, x=as.numeric(all.countries$west), y=as.numeric(all.countries$north), 
	LatLongDatum=rep("WGS84", NROW(all.countries)), Elevation=rep(0, NROW(all.countries)))

	addSite(Code="ANT", Name="Netherlands Antilles", x=0, y=0, 
	LatLongDatum="WGS84", Elevation=0)
	addSite(Code="YUG", Name="Yugoslavia", x=0, y=0, 
	LatLongDatum="WGS84", Elevation=0)
	addSite(Code="WORLD", Name="The World", x=0, y=0, LatLongDatum="WGS84", 
	Elevation=0)


	#############
	## Data Sources
	cat("Adding other Metadata\n")

	addISOMetadata(TopicCategory = "farming", Title = "FAOSTAT", Abstract = "FAOSTAT time-series and cross sectional data relating to food and agriculture. Detailed country-level data 
		on food consumption and production", ProfileVersion = "Unknown", MetadataLink = NULL)
	addSource(Organization="FAOSTAT", SourceDescription="FAOSTAT time-series and cross sectional data relating to food and agriculture", 
			SourceLink="faostat.fao.org", ContactName="LongExample Lissner",Metadata="FAOSTAT")
			
	addISOMetadata(TopicCategory = "society", Title = "IDP", Abstract = "Institutional Profiles Database: indicators on institutional characteristics of 123 countries", 
			ProfileVersion = "Unknown", MetadataLink = NULL)
		addSource(Organization="IPD", SourceDescription= "IDP database", SourceLink="www.cepii.fr/anglaisgraph/bdd/institutions.htm", ContactName="LongExample Lissner", Metadata = "IDP")
			
	addISOMetadata(TopicCategory = "society", Title = "IIASA", Abstract = "IIASA downscaled population projections for IPCC SRES", 
			ProfileVersion = "Unknown", MetadataLink = NULL)
		addSource(Organization="IIASA", SourceDescription= "IIASA downscaled population projections for IPCC SRES", SourceLink="http://ciesin.columbia.edu/datasets/downscaled/",
		ContactName="LongExample Lissner", Metadata = "IIASA")
			
	addISOMetadata(TopicCategory = "environment", Title = "WATCHproject", Abstract = "Results from the WATCH project", 
			ProfileVersion = "Unknown", MetadataLink = NULL)
		addSource(Organization="WATCH Project", SourceDescription="Results from the WATCH project", SourceLink="ftp://watch-r:wWread77@ftp.iiasa.ac.at/", ContactName="LongExample Lissner",Metadata="WATCHproject")

	addISOMetadata(TopicCategory = "health", Title = "WHO", Abstract = "Data on global health", ProfileVersion = "Unknown", MetadataLink = NULL)
		addSource(Organization="WHO", SourceDescription="Data on global health", SourceLink="http://apps.who.int/ghodata/#", ContactName="LongExample",Metadata="WHO")
		
	addISOMetadata(TopicCategory = "health", Title = "UNDP", Abstract = "Data on human development (Human Development Report)", 
			ProfileVersion = "Unknown", MetadataLink = NULL)
		addSource(Organization="UNDP", SourceDescription="Data on human development (Human Development Report)", SourceLink="hdrstats.undp.org/en", ContactName="LongExample",Metadata="UNDP")

	addISOMetadata(TopicCategory = "health", Title = "UNICEF", Abstract = "Childinfo - Monitoring the situation of women and children", 
			ProfileVersion = "Unknown", MetadataLink = NULL)	
		addSource(Organization="UNICEF", SourceDescription="Childinfo - Monitoring the situation of women and children", 
		SourceLink="http://www.childinfo.org/water_data.php", ContactName="LongExample", Metadata="UNICEF")
		
	addISOMetadata(TopicCategory = "utilitiesCommunication", Title = "IEA", Abstract = "International Energy Agency - Information on electricity and energy generation and use", 
			ProfileVersion = "Unknown", MetadataLink = NULL)		
		addSource(Organization="IEA", SourceDescription="International Energy Agency", SourceLink="http://www.iea.org/weo/electricity.asp", ContactName="LongExample", Metadata="IEA")



	#################
	### Quality Control Levels

	addQualityControlLevel(ID=6,Code="ok", Definition="The default")


	#############
	## Units


	addUnits(Name="cubic kilometers per year",Type="water availability",Abbreviation="km^3/yr")
	addUnits(Name="Kilocalories",Type="Energy",Abbreviation="kcal")
	addUnits(Name="m3/cap/yr",Type="Water availability",Abbreviation="m3/cap/yr")


	cat("Adding Synonyms\n")
	importSynonyms(system.file("longexample/syn.txt", package="RObsDat"))

	cat("Metadata prepared\n##############################\n")

	#############
	## Data
	#############

	###################
	# Total annual surface and subsurface run-off from WaterGAP; 1981-2010 and 2011-2040, Scenario A2, three models --> still in km^3/yr, needs to be calculated per cap

	cat("Importing water availability\n")
	wateravail.2000 <- read.csv(system.file("longexample/ofile_watergap_wfdnat_qstot_mean_1971_2000.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2030_ipsl <- read.csv(system.file("longexample/ofile_watergap_ipsla2nat_qstot_mean_2011_2040.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2030_echam <- read.csv(system.file("longexample/ofile_watergap_echama2nat_qstot_mean_2011_2040.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2030_cncm3 <- read.csv(system.file("longexample/ofile_watergap_cncm3a2nat_qstot_mean_2011_2040.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2060_ipsl <- read.csv(system.file("longexample/ofile_watergap_ipsla2nat_qstot_mean_2041_2070.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2060_echam <- read.csv(system.file("longexample/ofile_watergap_echama2nat_qstot_mean_2041_2070.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2060_cncm3 <- read.csv(system.file("longexample/ofile_watergap_cncm3a2nat_qstot_mean_2041_2070.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2090_ipsl <- read.csv(system.file("longexample/ofile_watergap_ipsla2nat_qstot_mean_2071_2100.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2090_echam <- read.csv(system.file("longexample/ofile_watergap_echama2nat_qstot_mean_2071_2100.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")
	wateravail.2090_cncm3 <- read.csv(system.file("longexample/ofile_watergap_cncm3a2nat_qstot_mean_2071_2100.csv", package="RObsDat"),header=TRUE,sep=",",stringsAsFactors=FALSE,na.strings="NA",dec=".")

	water.dat <- cbind(wateravail.2000,wateravail.2030_ipsl[,2],wateravail.2060_ipsl[,2],wateravail.2090_ipsl[,2],wateravail.2030_echam[,2],wateravail.2060_echam[,2],
	wateravail.2090_echam[,2],wateravail.2030_cncm3[,2],wateravail.2060_cncm3[,2],wateravail.2090_cncm3[,2])

	water.date <- strptime(c("1990-01-01","2030-01-01","2060-01-01","2090-01-01","2030-01-01","2060-01-01","2090-01-01","2030-01-01","2060-01-01","2090-01-01"), "%Y-%m-%d", tz="UTC") 

	water.vars <- c("water.available.cur","water.available.ipsl","water.available.ipsl","water.available.ipsl","water.available.echam","water.available.echam","water.available.echam",
	"water.available.cncm3","water.available.cncm3","water.available.cncm3")

	addCV("VariableName", term="water.available", definition="Total annual renewable water resources per country in km^3")
	addVariable(Code="water.available.cur", Name="water.available", Unit="km^3/yr") 
	addVariable(Code="water.available.ipsl", Name="water.available", Unit="km^3/yr")
	addVariable(Code="water.available.echam", Name="water.available", Unit="km^3/yr")
	addVariable(Code="water.available.cncm3", Name="water.available", Unit="km^3/yr")

	#Build synonym table
	getID("Site", water.dat[,1])

	#Import
	addDataValues(Date=water.date, Value=water.dat[,2:11],  Site = water.dat[,1], Variable = water.vars, Source = "WATCH Project", QualityControlLevel = "ok" )

	###################
	# Population projections SRES scenario A2 

	cat("Importing Population\n")
	population.dat <- read.csv2(system.file("longexample/Natl_Pop_Proj_A2.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, na.strings=c("na","nodata","-"), dec=".", encoding="utf8")
	population.date <- strptime(c("1990-01-01","1995-01-01","2000-01-01","2005-01-01","2010-01-01","2015-01-01","2020-01-01","2025-01-01","2030-01-01","2035-01-01",
				"2040-01-01","2045-01-01","2050-01-01","2055-01-01","2060-01-01","2065-01-01","2070-01-01","2075-01-01","2080-01-01","2085-01-01","2090-01-01","2095-01-01",
				"2100-01-01"),"%Y-%m-%d", tz="UTC") 
	population.vars <- c("population.A2")

	addCV("VariableName", term="population", definition="Projected number of inhabitants per country according to SRES")
	addVariable(Code="population.A2", Name="population", Unit="count") 

	getID("Site", population.dat[,1])

	#Import

	addDataValues(Date=population.date, Value=population.dat[,2:24],  Site = population.dat[,1], Variable = population.vars, Source = "IIASA", QualityControlLevel = "6" )

	################
	# Grand total: calorie availability/cap/day (2009 )from FAOSTAT 


	cat("Importing Calories\n")
	calories.dat  <- read.csv2(system.file("longexample/calorie_availability_2009_new.csv", package="RObsDat"), header=TRUE, dec=".", stringsAsFactors=FALSE, sep=";")
	calories.date <- strptime("2009-01-01", "%Y-%m-%d")
	calories.vars <- c("calories.total")

	addCV("VariableName", term="calories.total", definition="Total available calories/cap/day")
	addVariable(Code="calories.total", Name="Total available calories/cap/day", Unit="kcal") 

	#Build synonym table
	getID("Site", calories.dat[,1])

	#Import
		addDataValues(Date=calories.date, Value=calories.dat[,2],  Site = calories.dat[,1], Variable = calories.vars, Source = "FAOSTAT", QualityControlLevel = 6 )


	################
	# MDER: Minimum dietary requirements/cap/day

	cat("Importing  MDER\n")
	mder.dat  <- read.csv2(system.file("longexample/MinimumDietaryEnergyRequirement.csv", package="RObsDat"), header=TRUE, dec=".", stringsAsFactors=FALSE, sep=";")
	mder.date <- strptime("2008-01-01", "%Y-%m-%d")
	mder.vars <- c("mder")

	addCV("VariableName", term="mder", definition="Minimum dietary requirement: calories/cap/day")
	bla <- getMetadata("Variable", Name="dietary")
	print(bla)
	addVariable(Code="mder", Name="Minimum dietary requirement: calories/cap/day", Unit="kcal") 

	#Build synonym table
	getID("Site", mder.dat[,2])

	#Import
	addDataValues(Date=mder.date, Value=mder.dat[,7],  Site = mder.dat[,2], Variable = mder.vars, Source = "FAOSTAT", QualityControlLevel = 6 )

	   
	#########################
	##  IDP database, selected indicators: Population.participation;public.freedom_civil.society;stability.of.political.system;
	##       Domestic.public.security;Functioning.of.justice.system;Social.inclusion;geographic.coverage.of.public.services;Institutional.solidarity;
	##		 Traditional.solidarity;Micro.lending;Existence/absence.of.labour.legislation.and.measures;weak.employment.contract.rigidity
	######## 

	cat("Importing  Social Variables\n")
	IDP.dat <- read.csv2(system.file("longexample/IPD_2009_selected_indicators1.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=",")
	IDP.date <- strptime("2009-01-01",  "%Y-%m-%d", tz="UTC") 
	IDP.vars <- c(NA,"population.participation","public.freedom_civil.society","stability.of.political.system","Domestic.public.security","Functioning.of.justice.system",
	"Social.inclusion","geographic.coverage.of.public.services","Institutional.solidarity","Traditional.solidarity","Micro.lending","labour.legislation",
	"employment.contract")


	addCV("VariableName", term="Population.participation", definition="Population Participation")
	addCV("VariableName", term="public.freedom_civil.society", definition="Public Freedom and Civil Society")
	addCV("VariableName", term="stability.of.political.system", definition="Statibility of Political System")
	addCV("VariableName", term="Domestic.public.security", definition="Domestic Public Security")
	addCV("VariableName", term="Functioning.of.justice.system", definition="Functioning of Justice System")
	addCV("VariableName", term="Social.inclusion", definition="Social Inclusion")
	addCV("VariableName", term="geographic.coverage.of.public.services", definition="Geographic Coverage of Public Services")
	addCV("VariableName", term="Institutional.solidarity", definition="Institutional Solidarity")
	addCV("VariableName", term="Traditional.solidarity", definition="Traditional Solidarity")
	addCV("VariableName", term="Micro.lending", definition="Micro Lending")
	addCV("VariableName", term="labour.legislation", definition="Labour Legislation")
	addCV("VariableName", term="employment.contract", definition="Weak Employment Contracts")
	addVariable(Code="population.participation", Name="Population.participation", Unit="categorical") 
	addVariable(Code="public.freedom_civil.society", Name="public.freedom_civil.society", Unit="categorical") 
	addVariable(Code="stability.of.political.system", Name="stability.of.political.system", Unit="categorical") 
	addVariable(Code="Domestic.public.security", Name="Domestic.public.security", Unit="categorical") 
	addVariable(Code="Functioning.of.justice.system", Name="Functioning.of.justice.system", Unit="categorical") 
	addVariable(Code="Social.inclusion", Name="Social.inclusion", Unit="categorical") 
	addVariable(Code="geographic.coverage.of.public.services", Name="geographic.coverage.of.public.services", Unit="categorical") 
	addVariable(Code="Institutional.solidarity", Name="Institutional.solidarity", Unit="categorical") 
	addVariable(Code="Traditional.solidarity", Name="Traditional.solidarity", Unit="categorical") 
	addVariable(Code="Micro.lending", Name="Micro.lending", Unit="categorical") 
	addVariable(Code="labour.legislation", Name="labour.legislation", Unit="categorical") 
	addVariable(Code="employment.contract", Name="employment.contract", Unit="categorical") 

	#IDPcode.dat <- merge(IDP.dat, all.countries, by.x = "countrycode", by.y = "isoAlpha3")

	#Build synonym table
	getID("Site", IDP.dat[,1])

	#Import
	addDataValues(Date=IDP.date, Value=IDP.dat[,2:13],  Site = IDP.dat[,1], Variable = IDP.vars[2:13], Source = "IDP database", QualityControlLevel = 6)

	###################
	# Life expectancy at birth as used in HDI
	# added to "database_22_09_2011.db" and checked by LongExample

	cat("Importing  Life expect\n")
	lifeexp.dat <- read.csv2(system.file("longexample/life_expectancy_HDI_2009.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=";", dec=".")
	lifeexp.date <- strptime("2009-01-01",  "%Y-%m-%d", tz="UTC") 
	lifeexp.vars <- c("lifeexp")

	addCV("VariableName", term="lifeexp", definition="Life expectancy at birth")
	addVariable(Code="lifeexp", Name="lifeexp", Unit="yr") 

	#Build synonym table
	getID("Site", lifeexp.dat[,1])

	#Import
	addDataValues(Date=lifeexp.date, Value=lifeexp.dat[,2],  Site = lifeexp.dat[,1], Variable = lifeexp.vars, Source = "UNDP", QualityControlLevel = 6 )


	###################
	# PM10 concentration, lastest values: 2007 - 2010

	cat("Importing  PM10\n")
	pm10.dat <- read.csv2(system.file("longexample/pm10_latest_2009_Worldbank.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=",", dec=".")
	pm10.date <- strptime("2009-01-01",  "%Y-%m-%d", tz="UTC")
	pm10.vars <- c("pm10")

	addCV("VariableName", term="pm10", definition="Mean yearly PM10 concentrations")
	addVariable(Code="pm10", Name="pm10", Unit="ug/m^3") 

	#Build synonym table
	getID("Site", pm10.dat[,1])

	#Import
	addDataValues(Date=pm10.date, Value=pm10.dat[,2],  Site = pm10.dat[,1], Variable = pm10.vars, Source = "WHO", QualityControlLevel = 6 )


	###################
	# Solid fuel use

	cat("Importing  Solid Fuel\n")
	solid.fuel.dat <- read.csv2(system.file("longexample/solid_fuel_use2010.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=",", dec=".")
	solid.fuel.date <- strptime("2010-01-01",  "%Y-%m-%d", tz="UTC")
	solid.fuel.vars <- c("solid.fuel")

	addCV("VariableName", term="solid.fuel", definition="Mean yearly PM10 concentrations")
	addVariable(Code="solid.fuel", Name="solid.fuel", Unit="%") 

	#Build synonym table
	getID("Site", solid.fuel.dat[,1])

	#Import
	addDataValues(Date=solid.fuel.date, Value=solid.fuel.dat[,2],  Site = solid.fuel.dat[,1], Variable = solid.fuel.vars, Source = "WHO", QualityControlLevel = 6 )



	########################
	##  Water Data from Aquastat for Australia Analysis


	cat("Importing  Sanitation\n")
	watsan.dat <- read.csv2(system.file("longexample/access_sanitation2010.csv", package="RObsDat"), header=TRUE, na.strings=c("na","nodata","-"), stringsAsFactors=TRUE, sep=",")
	watsan.date <- strptime("2010-01-01",  "%Y-%m-%d", tz="UTC") 
	watsan.vars <- c("wat.total.improved","wat.piped","wat.other.improved","wat.other.unimproved","san.total.improved","san.shared","san.other.unimproved","san.open.def")

	addCV("VariableName", term="wat.total.improved", definition="Access to improved water sources (total)")
	addCV("VariableName", term="wat.piped", definition="Water piped on premises")
	addCV("VariableName", term="wat.other.improved", definition="Access to other improved water sources")
	addCV("VariableName", term="wat.other.unimproved", definition="Unimproved water sources")
	addCV("VariableName", term="san.total.improved", definition="Access to improved sanitation facilities")
	addCV("VariableName", term="san.shared", definition="Ahred sanitation facilities")
	addCV("VariableName", term="san.other.unimproved", definition="Unimproved sanitation access")
	addCV("VariableName", term="san.open.def", definition="Open defecation")
	addVariable(Code="wat.total.improved", Name="wat.total.improved", Unit="%") 
	addVariable(Code="wat.piped", Name="wat.piped", Unit="%") 
	addVariable(Code="wat.other.improved", Name="wat.other.improved", Unit="%") 
	addVariable(Code="wat.other.unimproved", Name="wat.other.unimproved", Unit="%") 
	addVariable(Code="san.total.improved", Name="san.total.improved", Unit="%") 
	addVariable(Code="san.shared", Name="san.shared", Unit="%") 
	addVariable(Code="san.other.unimproved", Name="san.other.unimproved", Unit="%") 
	addVariable(Code="san.open.def", Name="san.open.def", Unit="%") 

	getID("Site", watsan.dat[,1])
	#Build synonym table

	#Import
	addDataValues(Date=watsan.date, Value=watsan.dat[,2:9],  Site = watsan.dat[,1], Variable = watsan.vars, Source = "UNICEF", QualityControlLevel = 6 )


	#######################
	# Health care workforce --> is per 100000! needs to be converted

	cat("Importing  Health\n")

	HCW.dat <- read.csv(system.file("longexample/healthcare_workforce_sum1000_latest.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=",")
	HCW.date <- strptime(paste(HCW.dat[,2],"-01-01", sep=""),  "%Y-%m-%d", tz="UTC") 
	HCW.vars <- c("health.workers")

	addCV("VariableName", term="health.workers", definition="Density of health care work force (per 100000 population)")
	addVariable(Code="health.workers", Name="health.workers", Unit="#") 

	#Build synonym table
	getID("Site", HCW.dat[,1])

	#Import
	addDataValues(Date=HCW.date, Value=HCW.dat[,3],  Site = HCW.dat[,1], Variable = HCW.vars, Source = "WHO", QualityControlLevel = 6 )


	###################
	# Under 5 mortality
	# added to "database_22_09_2011.db" and checked by LongExample

	u5.dat <- read.csv2(system.file("longexample/under5mortality_rate.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=",")
	u5.date <- strptime("2010-01-01","%Y-%m-%d", tz="UTC") 
	u5.vars <- c("under5.mort")

	addCV("VariableName", term="under5.mort", definition="Probability of dying by age 5 per 1000 live births")
	addVariable(Code="under5.mort", Name="under5.mort", Unit="#") 


	#Build synonym table
	getID("Site", u5.dat[,1])

	#Import
		addDataValues(Date=u5.date, Value=u5.dat[,2],  Site = u5.dat[,1], Variable = u5.vars, Source = "WHO", QualityControlLevel = 6 )


	###################
	# Electrification rate
	# new version of table: countries with HDI very high set to 99, HDI high to 95
	# fully added to "database_22_09_2011.db" and checked by LongExample --> values updated by Lisei, added 2/12/2011

	elect.dat <- read.csv2(system.file("longexample/electricity_rate_2008_new.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=";", dec=".")
	elect.dat <- read.csv2(system.file("longexample/electrification_2009.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=",", dec=".")
	elect.date <- strptime("2008-01-01",  "%Y-%m-%d", tz="UTC") 
	elect.date <- strptime("2009-01-01",  "%Y-%m-%d", tz="UTC") 
	elect.vars <- c("electrate.total")

	addCV("VariableName", term="electrate.total", definition="Electrification rate of total population (percent of population with access to electricity)")
	addVariable(Code="electrate.total", Name="electrate.total", Unit="%") 


	#Build synonym table
	getID("Site", elect.dat[,1])

	#Import
		addDataValues(Date=elect.date, Value=elect.dat[,2],  Site = elect.dat[,1], Variable = elect.vars, Source = "IEA", QualityControlLevel = 6 )


	#########################
	##  Mean years of schooling and expected years of schooling as used in the new HDI
	######## fully added to "database_22_09_2011.db" and checked by LongExample


	schooling.dat <- read.csv2(system.file("longexample/years_of_schooling_HDI_2009.csv", package="RObsDat"), header=TRUE, stringsAsFactors=FALSE, sep=";", dec=".", na.strings="..")
	schooling.date <- strptime("2009-01-01",  "%Y-%m-%d", tz="UTC") 
	schooling.vars <- c("mean_schooling", "expected_schooling")

	addCV("VariableName", term="mean_schooling", definition="Mean Years of Schooling")
	addCV("VariableName", term="expected_schooling", definition="Expected mean Years of Schooling")
	addVariable(Code="mean_schooling", Name="mean_schooling", Unit="yr") 
	addVariable(Code="expected_schooling", Name="expected_schooling", Unit="yr") 

	 
	#Build synonym table
	getID("Site", schooling.dat[,1])


	#Import
	addDataValues(Date=schooling.date, Value=schooling.dat[,2:3],  Site = schooling.dat[,1], Variable = schooling.vars, Source = "UNDP", QualityControlLevel = 6 )


	####################

	req.vars.rest <- c("population.participation","public.freedom_civil.society","stability.of.political.system","Domestic.public.security","Functioning.of.justice.system",
					"Social.inclusion","Institutional.solidarity","Traditional.solidarity","Micro.lending","labour.legislation","employment.contract","calories.total",
					"mder","lifeexp","pm10","solid.fuel","wat.piped","wat.other.improved","wat.other.unimproved","san.total.improved",
					"san.shared","san.other.unimproved","san.open.def","health.workers","under5.mort","electrate.total","mean_schooling","expected_schooling")

	req.vars.water <- c("population.A2","water.available.cur","water.available.ipsl","water.available.echam","water.available.cncm3")
	
	all.water.dat <- getDataValues(Variable=req.vars.water)

#	all.data <- merge(all.water.dat, rest.data, by.x = "contryCode", by.y = "countryCode")

}
