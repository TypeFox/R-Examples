## ----results='hide', message=FALSE, warning=FALSE------------------------
# load the package
library(emuR)

create_emuDB(name = 'fromScratchDB', 
             targetDir = tempdir(), 
             verbose = F)

## ------------------------------------------------------------------------
# generate path to the empty fromScratchDB created above
dbPath = file.path(tempdir(), 'fromScratchDB_emuDB')
# load database
dbHandle = load_emuDB(dbPath, verbose = F)
print(dbHandle)

## ------------------------------------------------------------------------
add_levelDefinition(dbHandle, 
                    name = 'Phonetic', 
                    type = 'SEGMENT')

## ------------------------------------------------------------------------
list_levelDefinitions(dbHandle)

## ------------------------------------------------------------------------
summary(dbHandle)

## ------------------------------------------------------------------------
# add
add_levelDefinition(dbHandle,
                    name = 'Word',
                    type = 'ITEM')
# list
list_levelDefinitions(dbHandle)

## ------------------------------------------------------------------------
# list
list_attributeDefinitions(dbHandle, 
                          levelName = 'Phonetic')

## ------------------------------------------------------------------------
# add
add_attributeDefinition(dbHandle,
                        levelName = 'Phonetic',
                        name = 'SAMPA')
# list
list_attributeDefinitions(dbHandle,
                          levelName = 'Phonetic')

## ------------------------------------------------------------------------

ipaVowels = c('i', 'iː', 'u', 'uː', 'ə')

sampaVowels = c('i', 'i:', 'u', 'u:', '@')

# set legalLabels values for phonetic attributeDefinition
set_legalLabels(dbHandle,
                levelName = 'Phonetic',
                attributeDefinitionName = 'Phonetic',
                legalLabels = ipaVowels)

# get
get_legalLabels(dbHandle, 
                levelName = 'Phonetic', 
                attributeDefinitionName = 'Phonetic')

# set legalLabels values for phonetic attributeDefinition
set_legalLabels(dbHandle,
                levelName = 'Phonetic',
                attributeDefinitionName = 'SAMPA',
                legalLabels = sampaVowels)

# get
get_legalLabels(dbHandle, 
                levelName = 'Phonetic', 
                attributeDefinitionName = 'SAMPA')


## ------------------------------------------------------------------------
# add long vowels
add_attrDefLabelGroup(dbHandle,
                      levelName = 'Phonetic', 
                      attributeDefinitionName = 'Phonetic',
                      labelGroupName = 'long',
                      labelGroupValues = c('iː', 'uː'))

# add short vowels
add_attrDefLabelGroup(dbHandle,
                      levelName = 'Phonetic', 
                      attributeDefinitionName = 'Phonetic',
                      labelGroupName = 'short',
                      labelGroupValues = c('i', 'u', 'ə'))


# list
list_attrDefLabelGroups(dbHandle,
                        levelName = 'Phonetic',
                        attributeDefinitionName = 'Phonetic')


## ------------------------------------------------------------------------
# add
add_linkDefinition(dbHandle,
                   type = 'ONE_TO_MANY',
                   superlevelName = 'Word',
                   sublevelName = 'Phonetic')

## ------------------------------------------------------------------------
# get path to folder containing wav files 
# (in this case wav files that come with the wrassp package)
fp = system.file('extdata', package='wrassp')
# import media files into emuDB session called filesFromWrassp
import_mediaFiles(dbHandle, 
                  dir = fp, 
                  targetSessionName = 'filesFromWrassp', 
                  verbose = F)
# list session
list_sessions(dbHandle)
# list bundles
list_bundles(dbHandle)

## ------------------------------------------------------------------------
# show head of list_files
head(list_files(dbHandle))

## ---- results = "hide"---------------------------------------------------
# list all wav files in new emuDB
wavFilePaths = list.files(dbPath, 
                          pattern = "wav$", 
                          full.names = T, 
                          recursive = T)
# create folder to store zcr values in
outDirPath = file.path(tempdir(), 'zcranaVals')
dir.create(outDirPath)
# calculate zero-crossing-rate files
# using zcrana function of wrassp package
library(wrassp)
zcrana(listOfFiles = wavFilePaths, 
       outputDirectory = outDirPath)
# add zcr files to emuDB
add_files(dbHandle,
          dir = outDirPath, 
          fileExtension = 'zcr',
          targetSessionName = 'filesFromWrassp')

## ------------------------------------------------------------------------
# show head of list_files to check if files were added
head(list_files(dbHandle))

## ----results='hide', message=FALSE, warning=FALSE------------------------
# add track and calculate SSFF files by specifying 
# one of the signal processing functions the wrassp package provides
# (in this case the forest (formant estimation) function)
add_ssffTrackDefinition(dbHandle, 
                        name = 'formantValues', 
                        columnName = 'fm', 
                        fileExtension = 'fms',
                        onTheFlyFunctionName = 'forest')
# list
list_ssffTrackDefinitions(dbHandle)
# show head of list_files to check if files where added
head(list_files(dbHandle))

## ------------------------------------------------------------------------
# add 
add_ssffTrackDefinition(dbHandle, 
                        name = 'zeroCrossing', 
                        columnName = 'zcr', 
                        fileExtension = 'zcr')
# list
list_ssffTrackDefinitions(dbHandle)

## ------------------------------------------------------------------------
# list 
list_perspectives(dbHandle)

## ------------------------------------------------------------------------
# get order array of levels of default perspective
get_levelCanvasesOrder(dbHandle, 
                       perspectiveName = 'default')
# set order array of levels of default perspective
set_levelCanvasesOrder(dbHandle,
                       perspectiveName = 'default',
                       order = c('Phonetic'))

# get order array of levels of default perspective
get_levelCanvasesOrder(dbHandle, 
                       perspectiveName = 'default')


## ------------------------------------------------------------------------
# get order array of signals of default perspective
get_signalCanvasesOrder(dbHandle, 
                        perspectiveName = 'default')
# set order array of signals of default perspective
set_signalCanvasesOrder(dbHandle, 
                        perspectiveName = 'default',
                        order = c("OSCI", "SPEC", "formantValues", "zeroCrossing"))
# get order array of signals of default perspective
get_signalCanvasesOrder(dbHandle, 
                        perspectiveName = 'default')


## ----results='hide', message=FALSE, warning=FALSE------------------------
# create demo data in folder provided by the tempdir() function
create_emuRdemoData(dir = tempdir())
# get the path to a folder containing .wav & .TextGrid files that is part of the demo data
path2folder = file.path(tempdir(), "emuR_demoData", "TextGrid_collection")
# convert this TextGridCollection to the emuDB format
convert_TextGridCollection(path2folder, dbName = "myTGcolDB", 
                           targetDir = tempdir(), verbose = F)
# load database
dbHandle = load_emuDB(file.path(tempdir(), "myTGcolDB_emuDB"), verbose = F)

## ------------------------------------------------------------------------
# list levels
list_levelDefinitions(dbHandle)
# list ssffTracks
list_ssffTrackDefinitions(dbHandle)

## ------------------------------------------------------------------------
# add linkDefinition
add_linkDefinition(dbHandle, type = "ONE_TO_MANY",
                  superlevelName = "Syllable",
                  sublevelName = "Phoneme")

# list
list_linkDefinitions(dbHandle)

# envoke autobuild function
autobuild_linkFromTimes(dbHandle,
                       superlevelName = "Syllable",
                       sublevelName = "Phoneme",
                       convertSuperlevel = TRUE)

# list
list_levelDefinitions(dbHandle)

## ------------------------------------------------------------------------
# add linkDefinition
add_linkDefinition(dbHandle, type = "MANY_TO_MANY",
                  superlevelName = "Phoneme",
                  sublevelName = "Phonetic")

# list
list_linkDefinitions(dbHandle)

# invoke autobuild function
autobuild_linkFromTimes(dbHandle,
                       superlevelName = "Phoneme",
                       sublevelName = "Phonetic",
                       convertSuperlevel = TRUE)

# list
list_levelDefinitions(dbHandle)


## ----echo=FALSE, results='hide', message=FALSE, warning=FALSE------------
# disconnect to avoid file locking to sqliteDB that causes unlink
# to fail under windows
DBI::dbDisconnect(dbHandle$connection)

## ------------------------------------------------------------------------
# remove the newly generated emuDB and emuR_demoData as we will not be needing it 
# throughout the rest of this vignette
unlink(file.path(tempdir(), "emuR_demoData"), recursive = TRUE)
unlink(file.path(tempdir(),'myTGcolDB_emuDB'), recursive = TRUE)

