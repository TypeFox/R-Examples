
# ---------------------------------------------------------------------------
context("testing convert_BPFCollection")
# ---------------------------------------------------------------------------

sourceDirMain = file.path(tempdir(), "emuR_demoData")
testDir = file.path(tempdir(), "emuR_testthat")
dbName = "bpf_converter_test"

# Cleaning up (just in case)
unlink(file.path(testDir, dbName), recursive = T)

# ---------------------------------------------------------------------------
# Testing with original BPFs
# ---------------------------------------------------------------------------

sourceDir = file.path(sourceDirMain, "BPF_collection")
newDbFolderName = paste0(dbName, emuDB.suffix)
newDbPath = file.path(testDir, newDbFolderName) 
configPath = file.path(newDbPath, paste0(dbName, '_DBconfig.json')) 
            
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Code throws error when new levels are declared incorrectly",
          {
            # length(newLevels) != length(newLevelClasses)
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, newLevels = c("ABC"), newLevelClasses = c(1,2)),
                         regexp = "newLevels", ignore.case = T)
            
            # new level classes outside of range 1-5
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir,  dbName = dbName, verbose = F, newLevels = c("ABC"), newLevelClasses = c(6)),
                         regexp = "1.*5", ignore.case = T)
            
            # trying to change the class of an already existing BPF standard level
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, newLevels = c("ORT"), newLevelClasses = c(2)),
                         regexp = "standard", ignore.case = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Code throws error for failed directory checks",
          {
            # there is already a database of with the same name in the target dir
            dir.create(file.path(testDir, "something_silly_emuDB"))
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = "something_silly", verbose = F),
                         regexp = "directory.*already exists", ignore.case = T)
            unlink(file.path(testDir, "something_silly_emuDB"), recursive = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Error when using unifyLevels incorrectly.",
          {
            # unifyLevels without refLevel
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, unifyLevels = c("KAN")),
                         regexp = "unify.*reference", ignore.case = T)
            
            # refLevel in unifyLevels
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F,  refLevel = "ORT", unifyLevels = c("ORT", "KAN")),
                         regexp = "reference level", ignore.case = T)
            
            # class 2-5 level in unifyLevels
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, refLevel = "ORT", unifyLevels = c("GES")),
                         regexp = "unif", ignore.case = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Error when using refLevel incorrectly.",
          {
            # link-less refLevel
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, refLevel = "GES"),
                         regexp = "link-less.*reference level", ignore.case = T)
            
            # extractLevels on, but refLevel not in extractLevels
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, extractLevels = c("MAU", "TRN"),  refLevel = "ORT"),
                         regexp = "reference level", ignore.case = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Error when trying declare an unknown level in refLevel, extractLevels or unifyLevels.",
          {
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, extractLevels = c("ABC")),
                         regexp = "unknown level.*ABC", ignore.case = T)
            
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, refLevel = "ABC"),
                         regexp = "unknown level.*ABC", ignore.case = T)
            
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, refLevel = "ORT", unifyLevels = c("ABC")),
                         regexp = "unknown level.*ABC", ignore.case = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Error when segmentToEventLevels is used with a non-segment level",
          {
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, segmentToEventLevels = c("PRB")),
                         regexp = "segment", ignore.case = T)
            }
          )

# Cleaning up (just in case)
unlink(newDbPath, recursive = T)

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Conversion without reference level.",
          {
            convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F)
            
            # Format of data base.
            expect_true(newDbFolderName %in% list.dirs(testDir, full.names = F, recursive = F))
            expect_equal(length(list.files(newDbPath, recursive = F)), 2)
            expect_equal(length(list.files(file.path(newDbPath, "0000_ses"), recursive = F)), 7)
            expect_equal(length(list.files(file.path(newDbPath, "0000_ses", "msajc003_bndl"), recursive = F)), 2)
            
            # Correctness of config file.
            dbConfigLines = readLines(configPath, warn=F)
            dbConfig = jsonlite::fromJSON(paste(dbConfigLines, collapse=''), simplifyVector=F)
            
            # General & webAppConfig
            expect_equal(dbConfig$name, dbName)
            expect_equal(length(dbConfig$ssffTrackDefinitions), 0)
            expect_true(dbConfig$EMUwebAppConfig$activeButtons$saveBundle)
            expect_true(dbConfig$EMUwebAppConfig$activeButtons$showHierarchy)
            
            # Check that level canvas order is by order of appearance in BPF
            expect_equal(dbConfig$EMUwebAppConfig$perspectives[[1]]$levelCanvases$order, list("TRN", "MAU"))
                                    
            # Check that there are five level definitions (Utterance, KAN, ORT, TRN, MAU)
            expect_equal(length(dbConfig$levelDefinitions), 5)
                                    
            # Check that level names and types are correct
            expect_equal(sapply(dbConfig$levelDefinitions, function(x) x$name), c("Utterance", "KAN", "ORT", "TRN", "MAU"))
            expect_equal(sapply(dbConfig$levelDefinitions, function(x) x$type), c("ITEM", "ITEM", "ITEM", "SEGMENT", "SEGMENT"))
                                    
            # Check that each level has the appropriate amount of attribute definitions
            expect_equal(sapply(dbConfig$levelDefinitions, function(x) length(x$attributeDefinitions)), c(9, 1, 1, 1, 1))
            expect_equal(sapply(dbConfig$levelDefinitions, function(x) x$attributeDefinitions[[1]]$name), c("Utterance", "KAN", "ORT", "TRN", "MAU"))
                                    
            # Check that all header entries have become attributes of the Utterance level
            expect_equal(sapply(dbConfig$levelDefinitions[[1]]$attributeDefinitions, function(x) x$name),
                         c("Utterance", "LHD", "REP", "SNB", "SAM", "SBF", "SSB", "NCH", "SPN"))

            # No link definitions
            expect_equal(length(dbConfig$linkDefinitions), 0)
            
            # Correctness of one annot file (msajc003_annot)
            annotPath = file.path(newDbPath, "0000_ses", "msajc003_bndl", "msajc003_annot.json")
            dbAnnotLines = readLines(annotPath, warn=F)
            dbAnnot = jsonlite::fromJSON(paste(dbAnnotLines, collapse=''), simplifyVector=F)
                                    
            # Check that all levels have the appropriate number of items
            expect_equal(length(dbAnnot$levels[[1]]$items), 1)
            expect_equal(length(dbAnnot$levels[[2]]$items), 7)
            expect_equal(length(dbAnnot$levels[[3]]$items), 7)
            expect_equal(length(dbAnnot$levels[[4]]$items), 1)
            expect_equal(length(dbAnnot$levels[[5]]$items), 35)
                                    
            # Check individual items
            expect_equal(dbAnnot$levels[[1]]$items[[1]]$id, 1)
            expect_equal(dbAnnot$levels[[4]]$items[[1]]$sampleStart, 3800)
            expect_equal(dbAnnot$levels[[4]]$items[[1]]$sampleDur, 48199)
                                    
            # Check that all header entries have become labels of the Utterance item
            expect_equal(sapply(dbAnnot$levels[[1]]$items[[1]]$labels, function(x) x$name),
                         c("Utterance", "LHD", "REP", "SNB", "SAM", "SBF", "SSB", "NCH", "SPN"))
            expect_equal(sapply(dbAnnot$levels[[1]]$items[[1]]$labels, function(x) x$value),
                         c("msajc003", "Partitur 1.2.16", "unknown", "2", "20000", "01", "16", "1", "unknown"))
                                    
            # Check individual label
            expect_equal(dbAnnot$levels[[2]]$items[[3]]$labels[[1]]$value, "frendz")
                                    
            # Check that there are no links
            expect_equal(length(dbAnnot$links), 0)
            }
          )

# Cleaning up.
unlink(newDbPath, recursive = T)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Conversion with reference level.",
          {
            convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, refLevel = "ORT")
                 
            # Correctness of config file  
            dbConfigLines = readLines(configPath, warn=F)
            dbConfig = jsonlite::fromJSON(paste(dbConfigLines, collapse=''), simplifyVector=F)
            
            # Check that all link definitions are correct
            expect_equal(length(dbConfig$linkDefinitions), 5)
            expect_equal(sapply(dbConfig$linkDefinitions, function(x) x$superlevelName), c("ORT", "TRN", "ORT", "Utterance", "Utterance"))
            expect_equal(sapply(dbConfig$linkDefinitions, function(x) x$sublevelName), c("KAN", "ORT", "MAU", "ORT", "TRN"))
            expect_equal(sapply(dbConfig$linkDefinitions, function(x) x$type), c("ONE_TO_ONE", "ONE_TO_MANY", "ONE_TO_MANY", "ONE_TO_MANY", "ONE_TO_ONE"))
            
            # Correctness of one annot file (msajc003_bndl)
            annotPath = file.path(newDbPath, "0000_ses", "msajc003_bndl", "msajc003_annot.json")
            dbAnnotLines = readLines(annotPath, warn=F)
            dbAnnot = jsonlite::fromJSON(paste(dbAnnotLines, collapse=''), simplifyVector=F)
                                    
            # Check that Utterance item (ID 1) links to TRN item (ID 16) and ORT items (ID 9-15)
            expect_equal(unlist(sapply(dbAnnot$links, function(x) if(x$fromID == 1) x$toID)), c(9:16))
            
            # Check that TRN item (ID 16) links to ORT items (ID 9-15)
            expect_equal(unlist(sapply(dbAnnot$links, function(x) if(x$fromID == 16) x$toID)), c(9:15))
                                    
            # Check that ORT items (ID 9-15) link to KAN items (ID 2-8)
            expect_equal(unlist(sapply(dbAnnot$links, function(x) if(x$toID %in% c(2:8)) x$fromID)), c(9:15))
                                    
            # Check that ORT items (ID 9-15) link to MAU items (ID 17 and upwards)
            expect_equal(unique(unlist(sapply(dbAnnot$links, function(x) if(x$toID > 16) x$fromID))), c(9:15))
                                    
            # Check some individual links
            expect_equal(unlist(sapply(dbAnnot$links, function(x) if(x$fromID == 10) x$toID)), c(3,25))
            expect_equal(unlist(sapply(dbAnnot$links, function(x) if(x$fromID == 15) x$toID)), c(8, c(43:50)))
            }
          )

# Cleaning up
unlink(newDbPath, recursive = T)

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Conversion with unifyLevels",
          {
            convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, refLevel = "ORT", unifyLevels = c("KAN"))
            
            # Correctness of config file
            dbConfigLines = readLines(configPath, warn=F)
            dbConfig = jsonlite::fromJSON(paste(dbConfigLines, collapse=''), simplifyVector=F)
                                    
            # Check that there are only four levels (since KAN has become a label on ORT level)
            expect_equal(length(dbConfig$levelDefinitions), 4)
                                    
            # Check that ORT level has two labels ORT and KAN
            expect_equal(sapply(dbConfig$levelDefinitions[[2]]$attributeDefinitions, function(x) x$name), c("ORT", "KAN"))
                                    
            # Check that there is no link between ORT and KAN in link definitions
            expect_equal(sapply(dbConfig$linkDefinitions, function(x) x$superlevelName), c("TRN", "ORT", "Utterance", "Utterance"))
            expect_equal(sapply(dbConfig$linkDefinitions, function(x) x$sublevelName), c("ORT", "MAU", "ORT", "TRN"))
            expect_equal(sapply(dbConfig$linkDefinitions, function(x) x$type), c("ONE_TO_MANY", "ONE_TO_MANY", "ONE_TO_MANY", "ONE_TO_ONE"))
            
            # Correctness of one annot file (msajc003_bndl)
            annotPath = file.path(newDbPath, "0000_ses", "msajc003_bndl", "msajc003_annot.json")
            dbAnnotLines = readLines(annotPath, warn=F)
            dbAnnot = jsonlite::fromJSON(paste(dbAnnotLines, collapse=''), simplifyVector=F)
                                    
            # Check levels
            expect_equal(sapply(dbAnnot$levels, function(x) x$name), c("Utterance", "ORT", "TRN", "MAU"))
            expect_equal(dbAnnot$levels[[2]]$name, "ORT")
            expect_equal(dbAnnot$levels[[2]]$type, "ITEM")
                                    
            # Check that all items on level ORT have two labels, and that their names are ORT and KAN
            expect_equal(unique(sapply(dbAnnot$levels[[2]]$items, function(x) length(x$labels))), 2)
            expect_equal(unique(sapply(dbAnnot$levels[[2]]$items, function(x) x$labels[[1]]$name)), "ORT")
            expect_equal(unique(sapply(dbAnnot$levels[[2]]$items, function(x) x$labels[[2]]$name)), "KAN")
            
            # Check some individual labels
            expect_equal(dbAnnot$levels[[2]]$items[[3]]$labels[[1]]$value, "friends")
            expect_equal(dbAnnot$levels[[2]]$items[[3]]$labels[[2]]$value, "frendz")
            expect_equal(dbAnnot$levels[[2]]$items[[5]]$labels[[1]]$value, "was")
            expect_equal(dbAnnot$levels[[2]]$items[[5]]$labels[[2]]$value, "wQz")
            }
          )
# Cleaning up
unlink(newDbPath, recursive = T)

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Conversion with extractLevels.",
          {
            convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, extractLevels = c("MAU"))
            
            # Correctness of config file    
            dbConfigLines = readLines(configPath, warn=F)
            dbConfig = jsonlite::fromJSON(paste(dbConfigLines, collapse=''), simplifyVector=F)
                                    
            # Check that level definitions include only extractedLevels and Utterance
            expect_equal(length(dbConfig$levelDefinitions), 2)
            expect_equal(sapply(dbConfig$levelDefinitions, function(x) x$name), c("Utterance", "MAU"))
                                    
            # Check that there are no links defined (refLevel = NULL)
            expect_equal(length(dbConfig$linkDefinitions), 0)
           
            # Correctness of one annot file (msajc003_bndl)
            annotPath = file.path(newDbPath, "0000_ses", "msajc003_bndl", "msajc003_annot.json")
            dbAnnotLines = readLines(annotPath, warn=F)
            dbAnnot = jsonlite::fromJSON(paste(dbAnnotLines, collapse=''), simplifyVector=F)
                                    
            # Check that there are only two levels.
            expect_equal(sapply(dbAnnot$levels, function(x) x$name), c("Utterance", "MAU"))
                                    
            # Check that there are no links
            expect_equal(length(dbAnnot$links), 0)
            }
          )
# Cleaning up
unlink(newDbPath, recursive = T)

test_that("Loading emuDB",
          {
            convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, refLevel = "ORT", unifyLevels = c("KAN"))
            handle = load_emuDB(file.path(testDir, paste0(dbName, emuDB.suffix)), verbose = F)
            DBI::dbDisconnect(handle$connection)
            handle = NULL
          }
        )

# Cleaning up
unlink(newDbPath, recursive = T)


# ---------------------------------------------------------------------------
# Testing with manipulated BPFs
# ---------------------------------------------------------------------------

sourceDir = file.path(sourceDirMain, "BPF_collection_manipulated")
            
# Manipulated BPFs contain:
# msajc003.parmanipulated:
#       - multi-label label string on ORT level
#       - semicolon in link on KAN tier
#       - unknown level "XYZ" (with class 1 syntax)
#       - blank linkes
#       - missing SAM header tag
#       - segmental overlap on "MAU" tier
# msajc010.parmanipulated:
#       - empty BPF (no header or body)
# msajc012.parmanipulated:
#       - MAU -> ORT (ONE_TO_MANY)
            
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Correct call with necessary arguments",
          {
            convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, newLevels = c("XYZ"), 
                                           newLevelClasses = c(1), refLevel = "ORT", segmentToEventLevels = c("MAU"), bpfExt = "parmanipulated")
            
            # Correctness of config file
            dbConfigLines = readLines(configPath, warn=F)
            dbConfig = jsonlite::fromJSON(paste(dbConfigLines, collapse=''), simplifyVector=F)
                                    
            # Check that there are 6 levels defined (Utterance, KAN, ORT, TRN, MAU, XYZ)
            expect_equal(length(dbConfig$levelDefinitions), 6)
            
            # Check that MAU has been turned into an event level
            expect_equal(dbConfig$levelDefinitions[[5]]$type, "EVENT")
            expect_equal(dbConfig$levelDefinitions[[5]]$name, "MAU")
            
            # Check that new level XYZ has been added
            expect_equal(dbConfig$levelDefinitions[[6]]$name, "XYZ")
            expect_equal(dbConfig$levelDefinitions[[6]]$type, "ITEM")
            
            # Check that ORT has three label names defined (ORT, ABC, XYZ)
            expect_equal(sapply(dbConfig$levelDefinitions[[3]]$attributeDefinitions, function(x) x$name), c("ORT", "ABC", "XYZ"))
                                     
            # Check that link from ORT to MAU is MANY_TO_MANY
            expect_equal(dbConfig$linkDefinitions[[3]]$superlevelName, "ORT")
            expect_equal(dbConfig$linkDefinitions[[3]]$sublevelName, "MAU")
            expect_equal(dbConfig$linkDefinitions[[3]]$type, "MANY_TO_MANY")

            # Correctness of annot file msajc003_bndl
            annotPath = file.path(newDbPath, "0000_ses", "msajc003_bndl", "msajc003_annot.json")
            dbAnnotLines = readLines(annotPath, warn=F)
            dbAnnot = jsonlite::fromJSON(paste(dbAnnotLines, collapse=''), simplifyVector=F)
                              
            # Check that all labels on level 'MAU' have _start/_end suffix
            expect_true(all(unlist(sapply(dbAnnot$levels[[5]]$items, function(x) if(stringr::str_detect(x$labels[[1]]$value, "_start") || stringr::str_detect(x$labels[[1]]$value, "_end")) TRUE))))
            expect_equal(dbAnnot$levels[[5]]$items[[3]]$labels[[1]]$value, "@_start")
            expect_equal(dbAnnot$levels[[5]]$items[[4]]$labels[[1]]$value, "@_end")
                                    
            # Check that labels on level 'ORT' are correct
            expect_equal(sapply(dbAnnot$levels[[3]]$items[[2]]$labels, function(x) x$name), c("ABC", "XYZ"))
            expect_equal(sapply(dbAnnot$levels[[3]]$items[[2]]$labels, function(x) x$value), c("ABC_label", "XYZ_label"))
            expect_equal(sapply(dbAnnot$levels[[3]]$items[[3]]$labels, function(x) x$name), c("ORT"))
            expect_equal(sapply(dbAnnot$levels[[3]]$items[[3]]$labels, function(x) x$value), c("friends"))
                                    
            # Check that the item on 'KAN' with the semicolon does not have an incoming link
            expect_true(all(unlist(sapply(dbAnnot$links, function(x) if(x$toID == 5) FALSE))))

            # Correctness of annot file msajc010_bndl
            annotPath = file.path(newDbPath, "0000_ses", "msajc010_bndl", "msajc010_annot.json")
            dbAnnotLines = readLines(annotPath, warn=F)
            dbAnnot = jsonlite::fromJSON(paste(dbAnnotLines, collapse=''), simplifyVector=F)
            
            # Check that there is only one item (the utterance).
            expect_equal(sapply(dbAnnot$levels, function(x) length(x$items) > 0), c(T,F,F,F,F,F))
            
            # Check that the utterance item has only one label (the utterance name).
            expect_equal(length(dbAnnot$levels[[1]]$items[[1]]$labels), 1)
            expect_equal(dbAnnot$levels[[1]]$items[[1]]$labels[[1]]$name, "Utterance")
            expect_equal(dbAnnot$levels[[1]]$items[[1]]$labels[[1]]$value, "msajc010")
            
            # Check that there are no links.
            expect_equal(length(dbAnnot$links), 0)
            
            # Correctness of annot file msajc012_bndl
            annotPath = file.path(newDbPath, "0000_ses", "msajc012_bndl", "msajc012_annot.json")
            dbAnnotLines = readLines(annotPath, warn=F)
            dbAnnot = jsonlite::fromJSON(paste(dbAnnotLines, collapse=''), simplifyVector=F)
                                    
            # Check that links between ORT and MAU have been turned around (should have been MAU->ORT after parsing but ORT->MAU after link disambiguation)
            expect_equal(unlist(sapply(dbAnnot$links, function(x) if(x$toID == 21) x$fromID)), c(10,11,12))
            
            # Check that this annot file does not contain any items on the XYZ level
            expect_equal(sapply(dbAnnot$levels, function(x) length(x$items) > 0), c(T,T,T,T,T,F))
            }
          )
# Cleaning up
unlink(newDbPath, recursive = T)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Warnings (semicolon) are displayed if verbose.",
          {
            expect_warning(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = T, refLevel = "ORT", 
                                                         newLevels = c("XYZ"), newLevelClasses = c(1), segmentToEventLevels = c("MAU"), bpfExt = "parmanipulated"),
                           regexp = "between.*';'", ignore.case = T)
            }
          )
# Cleaning up
unlink(newDbPath, recursive = T)

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Conversion without overlap resolution on BPF with overlap causes error.",
          {
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F,  
                                                        newLevels = c("XYZ"), newLevelClasses = c(1), bpfExt = "parmanipulated"),
                         regexp = "overlap", ignore.case = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Conversion with unknown level name in a BPF causes error.",
          {
            expect_error(convert_BPFCollection(sourceDir = sourceDir, targetDir = testDir, dbName = dbName, verbose = F, bpfExt = "parmanipulated"),
                         regexp = "unknown level", ignore.case = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
test_that("Conversion with a mismatch between level class and BPF line causes error.",
          {
            expect_error(convert_BPFCollection(sourceDir = sourceDir,  targetDir = testDir, dbName = dbName, verbose = F, refLevel = "ORT", 
                                                        newLevels = c("XYZ"), newLevelClasses = c(5), segmentToEventLevels = c("MAU"), bpfExt = "parmanipulated"),
                         regexp = "class", ignore.case = T)
            }
          )
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Final clean-up (just in case)
unlink(newDbPath, recursive = T)

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
 
