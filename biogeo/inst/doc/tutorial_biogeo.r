## biogeo package tutorial
require(biogeo)
# Download the worldclim 10 minute bioclimatic variables data in generic grids format (.bil extension)
# Unzip the file into the extdata folder of the biogeo package

# You may want to download a better quality map of the world than the built-in "world" dataset
# see: www.thematicmapping.org/downloads/world_borders.php
# It can be read into R using the following:
# world <- readShapePoly("TM_WORLD_BORDERS-0.3.shp",IDvar="NAME") # Read the shape file

################################################################################################################
# 1. Data formatting for compatibility with biogeo
################################################################################################################
# Introduce the hypothetical species dataset
data(dat) # Access the species dataset
head(dat)  # View first few lines of data
data(world) # Access the country boundaries

# Check the data structure
(ck <- checkdatastr(dat)) # Check data structure
names(dat)

# Add fields necessary for running other functions
dat2 <- dat[,names(dat)[names(dat)!='ID']] #Remove the ID column
dat2 <- addmainfields(dat2, species='Species')
head(dat2)

# Add fields necessary for running other functions while specifying the names of required columns already in the dataset
data(gbifdat) # An example dataset from GBIF
names(gbifdat)
dat3 <- keepmainfields(gbifdat, Species='species', x='decimallongitude', y='decimallatitude')
head(dat3)
# You can also specifiy additional columns that you wish to keep
dat3 <- keepmainfields(gbifdat, Species='species', x='decimallongitude', y='decimallatitude', others=c('gbifid'))
head(dat3)
# Rename fields of dataset if you want to keep entire original dataset
names(gbifdat)
dat4 <- renamefields(gbifdat, ID='gbifid', x="decimallongitude", y="decimallatitude", Species="species")
names(dat4)

# Viewing formatted data
dev.new(width=5,height=4,noRStudioGD = TRUE)
a<-pointsworld(world, dat2) # basic plotting
a<-pointsworld(world, dat2, ext = c(10, 40,-36, -15)) # set the extent
a<-pointsworld(world, dat2, ext='p') # Zoom in on points
spU <- dat2[dat2$Species=='Species U',]
a<-pointsworld(world, spU, ext='p' ) # select a particular species

################################################################################################################
# 2. Convert coordinates to decimal degrees and find coordinates for localities that have no coordinates
################################################################################################################
# Introduce the "places" dataset for illustrating coordinate conversion
#places <- read.csv("Placestest.csv",stringsAsFactors=F) # If you want to read the csv file it can be found in the extdata folder
data(places) # Load the data
# Parse coordinates into separate columns for degrees, minutes and seconds
cdat <- dmsparse(places,x='long',y='lat',id='id') # parse of coordinates

# Find records that have decimal degree coordinates
fd <- finddecimals(places,x='long',y='lat')
places[which(fd==1),] # View these records in original dataset
cdat[which(fd==1),] # View these records in parsed dataset
fdms <- which(fd==0) # Select only those that are not in decimal degree format
places[fdms,]

# Dealing with missing coordinates in the dataset or incorrectly specified cardinal direction
places2 <- places # Create a new dataset
places2$long[1] <- "" # Assign an empty string to the longitude column for the first record
places2$long[2] <- NA # Assign a missing value
places2$long[5] <- '23 25 S' # Change the letter from E to S
cdat2 <- dmsparse(places2,x='long',y='lat',id='id') # Parse of coordinates
fe <- which(cdat2$exclude==1) #Select those with errors
cdat2[fe,]
fm <- missingcoords(places2$long,places2$lat) # Find row numbers for missing coordinates
cdat3 <- dmsparse(places2[-fm,],x='long',y='lat',id='id') # Remove the rows with missing coords

# When there are no delimiters between degrees, minutes and seconds
coordstr <- '2344E'  # degrees, minutes and decimal seconds
fmt <- 'ddmmL'
dmsabs(coordstr,fmt)
# When there are only spaces between degrees, minutes and seconds
coordstr <- '23 44 25.3 E'  # degrees, minutes and decimal seconds
fmt <- 'dd mm ss.s L'
dmsabs(coordstr,fmt)
# Degrees and minutes separated by a fullstop
coordstr <- '44.25 E' 
fmt <- 'dd.mm L'
dmsabs(coordstr,fmt)
# Decimal degrees
coordstr <- '-25.19782' 
fmt <- 'dd.ddddd'
dmsabs(coordstr,fmt)

# Get the format of each of the coordinates
coords <- c('44 25 E','21.20 E','W14.03','12.35.16 E','09.26.08 W')
format <- getformat(coords) #  Get the format of each of the coordinates
data.frame(coords,format)

# Get the unique list of coordinates
coords <- c('44 25 E','21 20 E','W14.03','12.35.16 E','09.26.08 W')
(fmx<-uniqueformats(coords)) # List of unique formats

# Mixed formats in the same dataset
xtxt <- c('44 25 E','12.35.16 E','21.20 E','14.03 E','09.26.08 W') # x-coordinates
fmx <- uniqueformats(xtxt) # List of unique formats
a <- as.character(fmx[,2])
fmtstr <- c("dd*mm*L" ,"dd.mm.ss*L","dd.mm*L")
parsecoords(xtxt,fmx,fmtstr) # parse the x-coordinates

# Convert decimal degrees into degrees, minutes and seconds (for longitude)
(dmslong <- dd2dmslong(dat2$x[1:5])) # Decimal degrees of longitude to dms
(dmslat <- dd2dmslat(dat2$y[1:5])) # Decimal degrees of latitude to dms

# Convert dms to decimal degrees (for longitude)
dec <- dms2dd(cdat$xdeg,cdat$xmin,cdat$xsec,cdat$EW) 
data.frame(xdeg=cdat$xdeg,xmin=cdat$xmin,xsec=cdat$xsec,EW=cdat$EW,dec)

# Extract coordinates from Google Earth
data(dat) # Access the species dataset
data(world) # Access the country boundaries
pointsworld(world, dat, ext=c(17.5,20.5,-34.7,-34)) # Plot the problem point (zoomed in)
dat[dat$ID==689,] # View the coordinates and locality name for this point before updating
# Move pointer over desired location (Gordons Bay) in Google Earth and then CTRL-SHIFT-C to copy coords to clipboard
dat <- fromGEarth(dat,ID=689) # Get the coordinates from the clipboard
dat[dat$ID==689,] # view the coordinates for this point after updating the coords

################################################################################################################
# 3. Identify duplicate records to prevent pseudoreplication
################################################################################################################
data(dat)
dat2 <- duplicatesexclude(dat,res=20)
dat2[dat2$Reason=='Duplicated',] # View duplicated records

################################################################################################################
# 4. Identify records that may be too imprecise for the analysis
################################################################################################################
data(dat) # Access the species dataset
names(dat) # Column names of species dataset

#Check precision
datpc <- precisioncheck(dat, x='x', y='y', s=10, e=60)
datpc[datpc$preci==1,] #View records with possible precision problems

#Check precision in relation to environmental raster
data(dem)
datpce <- precisionenv(dat, dem, x='x', y='y')
datpce[datpce$envpreci==1,] #View records with possible precision problems

################################################################################################################
# 5. Identify records that likely have incorrect coordinates using geographical and environmental information.
################################################################################################################
data(dat)
data(datm)
world <- readShapePoly('TM_WORLD_BORDERS-0.3.shp',IDvar="NAME")
data(dem)


#---------------------------------------------------------------------------------------------------------------
# 5.1. Identify points that are in the wrong environment (e.g. terrestrial species in the sea OR marine species on land)
#---------------------------------------------------------------------------------------------------------------
# Move points in the sea that are close to land to the nearest cell on land (for terrestrial species)
# View points in the wrong environment using pointsworld
pointsworld(world, dat, ext="p") #Points in the wrong environment are shown in red, those in the correct environment in blue
datx <- nearestcell(dat,dem) # Returns a list with (1) the modified data and (2) the IDs and coordinates of moved points
s1 <- datx$moved # The points moved
dat2 <- datx$dat # The modified data
cr <- datx$dat$Correction
s2 <- which(str_detect(cr,"7"))
datx$dat[s2,] #View data for records that were moved
pointsworld(world, dat, ext="p") # Plots points on world map and checks for errors
points(s1$x,s1$y,pch=18) # Display the moved points in black

# Move points on land that are close to the sea to the nearest cell in the sea (for marine species)
# View points in the wrong environment using pointsworld
rstm<-raster(xmn=-180, xmx=180, ymn=-90, ymx=90,res=10/60,vals=1)#
data(msk10) # load indices of land cells
rstm[msk10]<-NA # assign NA to the land cells
pointsworld(world, datm, ext="p") #Points in the wrong environment are shown in red, those in the correct environment in blue
datmx <- nearestcell(datm, rstm) # Returns a list with (1) the modified data and (2) the IDs and coordinates of moved points
s1 <- datmx$moved # The points moved
datm2 <- datmx$dat #The modified data
cr <- datmx$dat$Correction
s2 <- which(str_detect(cr,"7"))
datmx$dat[s2,] #View data for records that were moved
pointsworld(world, datm,ext="p") # Plots points on world map and checks for errors
points(s1$x,s1$y,pch=18) # Display the moved points in black

# Identify points that are in the wrong environment (e.g. terrestrial species in the sea OR marine species on land)
d2 <- missingvalsexclude(dem, dat) # identify records with missing values
d2[d2$Exclude==1,] # Exclude = 1 are records with NA values of rst

#---------------------------------------------------------------------------------------------------------------
# 5.2. Get alternative coordinates using geographical visualization
#---------------------------------------------------------------------------------------------------------------
# View alternative coordinates for a selected record. For all species
d3 <- alternatives(dat,group1="Species",group2="",world,dem,locality="LocalityName",pos="bottomleft",ext="p")
# View alternative coordinates for a selected record. For one species
d4 <- alternatives2(dat,g1="Species E",group1="Species",group2="",world,dem,locality="LocalityName",
                    pos="bottomright", ext=c(-28,60,-40,40))

#---------------------------------------------------------------------------------------------------------------
# 5.3. Get alternative coordinates using environmental variables
#---------------------------------------------------------------------------------------------------------------
# First get environmental data into a raster stack 
fd <- system.file(package="biogeo") # Find path for biogeo package
foldenv <- file.path(fd,"extdata", fsep = .Platform$file.sep) # The ex folder
ev <- env2stack(foldenv, vars = c("bio1","bio12","bio5","bio6"), fext="bil") # Stack
env <- extract(ev, cbind(dat$x,dat$y))
edat <- cbind(dat,env)
plotsetup(6,6) # Open windows for plotting (try plotsetupx11 if you are not using Windows)
g1 = "Species U" # Select species
vars = c("bio1","bio12") #Select environmental variables
wclim(f=c(1,12)) #View names of environmental variables
d5 <- alternativesenv(edat,g1,group1="Species",ev,vars,world,xname="Annual Mean Temperature",
                      yname="Annual Precipitation",dem,locality="LocalityName",ext="p")

#---------------------------------------------------------------------------------------------------------------
# 5.4. Alternatives using either geographical or environmental information
#---------------------------------------------------------------------------------------------------------------
# Geographical and environmental space (two variables) - identify points
sp <- "Species U" # Choose species
plotsetup(6,6) # Open plotting windows (see also plotsetupx11())
d6 <- geo2envid(edat,sp,group1="Species",group2="",world,xc="bio1",yc="bio5",xname="Ann. Precip.",
                yname="Ann. Mean Temp.",showrecord="",ext="p")
f <- which(d6$ID%in%c(1981,1983)) # Show record(s) selected using geo2envid
d6[f,]

# Geographical and environmental space (principal components) - identify points
sp <- "Species U" # Choose species
plotsetup(6,6) # Open plotting windows (see also plotsetupx11())
envvars <- c("bio1","bio12","bio5","bio6") #Select environmental variables for PCA
d7 <- geo2envpca(edat,sp,group1="Species",group2="",world,scaling=1,vars=envvars,
               showrecord="1981",ext="p")
f <- which(d7$ID%in%c(1981)) # Show record(s) selected using geo2envid
d7[f,]

#---------------------------------------------------------------------------------------------------------------
# 5.5. errorcheck provides a suite of error detection tools
#---------------------------------------------------------------------------------------------------------------
d8 <- errorcheck(world, dem, dat=edat, countries="Country", countryfield="NAME", 
                 vars=c("bio1", "bio12", "bio5", "bio6"), res=10,elevc="elev",diff=50)
# Duplicated records (within grid cells of 'dem')
d8[d8$dups==1,]
# Records for which the recorded country does not match that based on the record's coordinates
d8[d8$CountryMismatch==1,]
# Records in the wrong environment (e.g. terrestrial species with records in the sea)
d8[d8$wrongEnv==1,]
# Records with low precision
d8[d8$lowprec==1,]
# Records with elevation mismatch
d8[d8$elevMismatch==1,]
# Records with errors for any of the above
d8[d8$error==1,]
# Records with environmental outliers
d8[d8$spperr==1,]

# An example using errorcheck to illustrate plotting of outliers
d8 <- errorcheck(world, dem, dat=edat, countries="Country", countryfield="NAME", 
                 vars=c("bio1", "bio2", "bio5", "bio6"), res=10,elevc="",diff=50)
sp <- "Species U"
plotsetup(6,6) # 
d9 <- geo2envid(d8, g1=sp, group1="Species", group2="bio6_j", world, xc="bio1", yc="bio5", xname="Ann. Precip.",
                yname="Ann. Mean Temp.", showrecord="", ext="p") #Outliers are plotted in blue


#---------------------------------------------------------------------------------------------------------------
# 5.6. Check that elevation (or depth) of record matches that of values extracted using coordinates
#---------------------------------------------------------------------------------------------------------------
data(gbifdat) # Get example data with already recorded altitudes
data(dem) # Get altitude raster
gb <- keepmainfields(gbifdat,ID='',Species='species',x='decimallongitude',y='decimallatitude',
                     others=c('gbifid','elevation')) # Convert example data to biogeo format
gb <- gb[ gb$Species=='Heterotheca villosa', ] # Keep data for only one species
gba<-elevcheck(gb,dem,elevc="elevation",diff=50)
f<-which(gba$elevMismatch==1)
gbifdat[f,] # View original data for records with altitude mismatches

# quick clean
data(msk10) # load a 10 minute mask index
dat2<-quickclean(world,dat,ID='ID',Species='Species',x='x',y='y',countries = "",others='',res=10,msk=msk10,ext="")


################################################################################################################
# Data summaries and output
################################################################################################################
data(dat)
data(dem)

# View names of bioclim variables
wclim(full=T) # This gives the names of all the bioclim variables from WorldClim
wclim(f=c(1,3,9)) # This gives the names for specific bioclim variables, i.e. bio1, bio3 and bio9

# To count the number of records per species in a dataset
dat2 <- duplicatesexclude(dat,res=10) # Get duplicate records
(nsp <- speciescount(dat2, orderby="nuniq")) # Number of records per species. Ordered by the number of unique records

# To view records that were modified during the cleaning process
datx <- nearestcell(dat, dem) # Move records in the sea to the nearest land cell (for terrestrial species)
dat3 <- datx$dat
f <- modifiedtoday(dat3)
dat3[f,] # Extract the records that were modified today
f2 <- modified(dat3,"01-01-2015 00:00:00","31-12-2015 00:00:00") # Records modified between two dates/times
dat3[f2,]

# Richness maps
data(dat)
data(dem)
dem2<-crop(dem,c(15,35,-36,-23))
rich <- richnessmap(dat, dem2, option="richness")
colPal <- colorRampPalette(colors=c('#556270','#4ECDC4','#C7F464','#FF6B6B','#C44D58'))
plot(rich, col=colPal(100))
nrec <- richnessmap(dat,dem,option="records")
writeRaster(rich, filename="richtest.asc", datatype='INT4S', overwrite=TRUE)
writeRaster(nrec, filename="nrectest.asc", datatype='INT4S', overwrite=TRUE)

# richness maps where spatial resolution can be specified without a raster
ex1 <- c(15,35,-36,-23)
rich<-richness(dat,res=20,option="richness",buf=5,ext=ex1)
colPal <- colorRampPalette(colors=c('#556270','#4ECDC4','#C7F464','#FF6B6B','#C44D58'))
plot(rich, col=colPal(100))

# quick richness map
data(dat)
data(msk20) # load a 20 minute mask index
ex1 <- c(15,35,-36,-23) # set the extent
rich<-quickrich(dat,ID='ID',Species='Species',x='x',y='y',countries = "",others='',res=20,msk=msk20,ext=ex1)
plot(rich, col=palette(gray(seq(0,.9,len = 25))))


# Exporting data
# Save to a kml that can be opened in Google Earth
ss <- sample(1:nrow(dat),size=50)
dat2 <- dat[ss,] # Randomly sample 50 records
xy <- data.frame(dat2$x,dat2$y)
kmlEx <- SpatialPointsDataFrame(xy,data=dat2,proj4string=CRS("+proj=longlat +datum=WGS84"))
require(plotKML)
icon = "http://maps.google.com/mapfiles/kml/pal2/icon18.png" # Select an icon
# The code below will create a kml file called "kmlEx.kml" in the working directory that you can open in Google Earth
kml(kmlEx, shape=icon, colour=Species, labels=ID, size=1) 

# Write data to a point shapefile (for use in a GIS package)
fn <- "shapetest.shp" # The name of the shapefile
points2shape(dat, x="x", y="y", fn=fn) # This file can now be opened in a GIS package

# Write to a comma separated values file (csv file) - can read using MS Excel
write.csv(dat, file="csvtest.csv", row.names=F)

################################################################################################################
# A practical example using data from GBIF
################################################################################################################
# Download some data from GBIF
library(rgbif) # Load the rgbif package
(keys <- name_suggest(q='Panthera leo', rank='species', fields=c('key','canonicalName'))) # Search for correct name
key <- keys$key[1] # Choose the first key
lionsGBIF <- occ_search(taxonKey=key, hasCoordinate=T, limit=10000, fields='all') # Download data from GBIF
lionsGBIF$hierarchy # Check that the correct species has been downloaded
lions <- lionsGBIF$data # Get only the data

# Some basic cleaning of GBIF data
names(lions)
fields <- c('name','key','basisOfRecord','scientificName','decimalLongitude','decimalLatitude','year','issues',
            'countryCode','country','gbifID','institutionCode','occurrenceRemarks','locality','elevation','sex',
            'coordinateAccuracy','fieldNotes')
lions2 <- lions[,names(lions)%in%fields] # Keep only wanted fields
world <- readShapePoly('TM_WORLD_BORDERS-0.3.shp')
windows(7,5)
pointsworld(world, dat=lions2, x='decimalLongitude', y='decimalLatitude') #View uncleaned records

# Remove fossil records
summary(as.factor(lions2$basisOfRecord))
lions2 <- lions2[lions2$basisOfRecord!='FOSSIL_SPECIMEN',]
pointsworld(world, dat=lions2, x='decimalLongitude', y='decimalLatitude') #View records without fossil specimens

# Sometimes you may want to check species names for synonomy issues
unique(lions2$scientificName) # For example, you can check these against the ITIS database: http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=183803

# Get data into biogeo format
bglions <- keepmainfields(lions2, ID='key', Species='name', x='decimalLongitude', y='decimalLatitude',
                          others=c('basisOfRecord','scientificName','year','issues','countryCode','country',
                                   'gbifID','institutionCode','occurrenceRemarks','locality','elevation','sex',
                                   'coordinateAccuracy','fieldNotes'))

# Get environmental data
fd <- system.file(package="biogeo") # Find path for biogeo package
foldenv <- file.path(fd,"extdata", fsep = .Platform$file.sep) # The ex folder
wclim(full=T)
ev <- env2stack(foldenv, vars = c("bio1","bio12","bio5","bio6"), fext="bil") # Stack
env <- extract(ev, cbind(bglions$x,bglions$y))
bglions <- cbind(bglions,env)

# Use errorcheck to perform multiple corrections
bglions <- errorcheck(world, dem, dat=bglions, countries="country", countryfield="NAME", 
                      vars=c("bio1", "bio12", "bio5", "bio6"), res=10)

# Identify points that are in the sea
bglions <- missingvalsexclude(dem, bglions) # Identify records with missing values
nrow(bglions[bglions$Reason=='No raster values',]) # Number of records that are in the sea
bglions[bglions$Reason=='No raster values',] # View records in the sea

# Move points in the sea that are close to land to the nearest cell on land
bglionsx <- nearestcell(bglions, dem) # None of the records are close enough to a land cell to justify moving them

# View records in geographical and environmental space and decide whether to exclude any records
plotsetup(6,6) # see also plotsetupx11()
envvars <- c("bio1","bio12","bio5","bio6") #Select environmental variables for PCA
bglions <- geo2envpca(bglions, g1='Panthera leo', group1="Species", group2="bio1_j", world, scaling=1, vars=envvars, ext="p")
# You will notice that all the European and some North American records cluster together and are also outliers