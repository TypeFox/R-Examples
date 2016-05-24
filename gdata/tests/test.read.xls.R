library(gdata)

if ( ! 'XLSX' %in% xlsFormats() )
  {
    try( installXLSXsupport() )
  }

# iris.xls is included in the gregmisc package for use as an example
xlsfile <- file.path(path.package('gdata'),'xls','iris.xls')

iris.1 <- read.xls(xlsfile) # defaults to csv format
iris.1

iris.2 <- read.xls(xlsfile,method="csv") # specify csv format
iris.2

iris.3 <- read.xls(xlsfile,method="tab") # specify tab format
iris.3

stopifnot(all.equal(iris.1, iris.2))
stopifnot(all.equal(iris.1, iris.3))

exampleFile <- file.path(path.package('gdata'),'xls',
                         'ExampleExcelFile.xls')

exampleFileX <- file.path(path.package('gdata'),'xls',
                         'ExampleExcelFile.xlsx')

# see the number and names of sheets:
sheetCount(exampleFile)

if(! 'XLSX' %in% xlsFormats() )
  {
    cat("************************************************************\n")
    cat("** DIFF IN THIS SECTION IS EXPECTED BECAUSE PERL PACKAGES **\n")
    cat("** FOR SUPPORTING XLSX ARE NOT INSTALLED                  **\n")
    cat("************************************************************\n")
  } else {
    sheetCount(exampleFileX)
  }


sheetNames(exampleFile)

if(! 'XLSX' %in% xlsFormats() )
  {
    cat("************************************************************\n")
    cat("** DIFF IN THIS SECTION IS EXPECTED BECAUSE PERL PACKAGES **\n")
    cat("** FOR SUPPORTING XLSX ARE NOT INSTALLED                  **\n")
    cat("************************************************************\n")
  } else {
    sheetNames(exampleFileX)
  }


example.1 <- read.xls(exampleFile, sheet=1) # default is first worksheet
example.1

example.2 <- read.xls(exampleFile, sheet=2) # second worksheet by number
example.2

example.3 <- read.xls(exampleFile, sheet=3, header=FALSE) # third worksheet by number
example.3

example.4 <- read.xls(exampleFile, sheet=4, header=FALSE) # fourth worksheet by number
example.4

if(! 'XLSX' %in% xlsFormats() )
  {
    cat("************************************************************\n")
    cat("** DIFF IN THIS SECTION IS EXPECTED BECAUSE PERL PACKAGES **\n")
    cat("** FOR SUPPORTING XLSX ARE NOT INSTALLED                  **\n")
    cat("************************************************************\n")
  } else {
    example.x.1 <- read.xls(exampleFileX, sheet=1) # default is first worksheet
    print(example.x.1)

    example.x.2 <- read.xls(exampleFileX, sheet=2) # second worksheet by number
    print(example.x.2)

    example.x.3 <- read.xls(exampleFileX, sheet=3, header=FALSE) # third worksheet by number
    print(example.x.3)

    example.x.4 <- read.xls(exampleFileX, sheet=4, header=FALSE) # fourth worksheet by number
    print(example.x.4)

    data <- read.xls(exampleFileX, sheet="Sheet Second") # and by name
    print(data)

    # load the third worksheet, skipping the first two non-data lines...
    data <- read.xls(exampleFileX, sheet="Sheet with initial text", skip=2)
    print(data)
  }

## Check handling of skip.blank.lines=FALSE

example.skip <- read.xls(exampleFile, sheet=2, blank.lines.skip=FALSE)
example.skip

if(! 'XLSX' %in% xlsFormats() )
  {
    cat("************************************************************\n")
    cat("** DIFF IN THIS SECTION IS EXPECTED BECAUSE PERL PACKAGES **\n")
    cat("** FOR SUPPORTING XLSX ARE NOT INSTALLED                  **\n")
    cat("************************************************************\n")
  } else {
    example.x.skip <- read.xls(exampleFileX, sheet=2, blank.lines.skip=FALSE)
    example.x.skip
  }



## Check handing of fileEncoding for latin-1 characters

latin1File  <- file.path(path.package('gdata'),'xls', 'latin-1.xls' )
latin1FileX <- file.path(path.package('gdata'),'xls', 'latin-1.xlsx')

if(.Platform$OS.type=="unix")
  {
      example.latin1 <- read.xls(latin1File,
                                 fileEncoding='latin1',
                                 encoding='latin1',
                                 stringsAsFactors=FALSE)
  } else {
      example.latin1 <- read.xls(latin1File,
                                        #fileEncoding='latin1',
                                 encoding='latin1',
                                 stringsAsFactors=FALSE)
  }

if(! 'XLSX' %in% xlsFormats() )
  {
    cat("************************************************************\n")
    cat("** DIFF IN THIS SECTION IS EXPECTED BECAUSE PERL PACKAGES **\n")
    cat("** FOR SUPPORTING XLSX ARE NOT INSTALLED                  **\n")
    cat("************************************************************\n")
  } else {
      if(.Platform$OS.type=="unix")
          {
              example.latin1.x <- read.xls(latin1FileX,
                                         fileEncoding='latin1',
                                         encoding='latin1',
                                         stringsAsFactors=FALSE)
          } else {
              example.latin1.x <- read.xls(latin1FileX,
                                        #fileEncoding='latin1',
                                         encoding='latin1',
                                         stringsAsFactors=FALSE)
          }
  }


## Check handling of very wide file

wideFile  <- file.path(path.package('gdata'),'xls', 'wide.xls' )
wideFileX <- file.path(path.package('gdata'),'xls', 'wide.xlsx')

example.wide <- read.xls(wideFile)
stopifnot(dim(example.wide)==c(0,256))

if( !'XLSX' %in% xlsFormats() )
  {
    cat("************************************************************\n")
    cat("** DIFF IN THIS SECTION IS EXPECTED BECAUSE PERL PACKAGES **\n")
    cat("** FOR SUPPORTING XLSX ARE NOT INSTALLED                  **\n")
    cat("************************************************************\n")
  } else {
    example.wide.x <- read.xls(wideFileX)
    stopifnot(dim(example.wide.x)==c(0,16384))
  }

## Check handling of files with dates calulcated relative to
## 1900-01-01 and 1904-01-01

file.1900  <- file.path(path.package('gdata'),'xls', 'ExampleExcelFile_1900.xls' )
file.1904  <- file.path(path.package('gdata'),'xls', 'ExampleExcelFile_1904.xls' )
fileX.1900 <- file.path(path.package('gdata'),'xls', 'ExampleExcelFile_1900.xlsx')
fileX.1904 <- file.path(path.package('gdata'),'xls', 'ExampleExcelFile_1904.xlsx')

example.1900 <- read.xls(file.1900, sheet=3, header=FALSE)
example.1900

example.1904 <- read.xls(file.1904, sheet=3, header=FALSE)
example.1904

exampleX.1900 <- read.xls(file.1900, sheet=3, header=FALSE)
exampleX.1900

exampleX.1904 <- read.xls(file.1904, sheet=3, header=FALSE)
exampleX.1904

# all colmns should be identical
stopifnot( na.omit(example.1900  == exampleX.1900) )
stopifnot( na.omit(example.1904  == exampleX.1904) )

# column 8 will differ by 1462 due to different date baselines (1900 vs 1904)
stopifnot( na.omit(example.1900 [,-8] == example.1904 [,-8]) )
stopifnot( na.omit(exampleX.1900[,-8] == exampleX.1904[,-8]) )

stopifnot( na.omit(example.1900 [,8] - example.1904 [,8]) == 1462 )
stopifnot( na.omit(exampleX.1900[,8] - exampleX.1904[,8]) == 1462 )
