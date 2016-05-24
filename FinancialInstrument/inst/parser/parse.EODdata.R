# Script for parsing downloaded Global Index data from
# http://EODdata.com (which requires purchase),
# and instructions for defining the symbols as instruments and 
# registering them with quantmod's getSymbols function.

# Peter Carl 

# This script requires the following directory structure:
# filesroot [directory set in the script below]
#   Each symbol's processed csv files are stored in sub-directories
#   named for each symbol, e.g., ~/Data/EOD\ Global\ Indexes/XMI.IDX.
#   These directories and files will be created and updated by this script.
# filesroot/.incoming
#   New or updated zip files should be placed here for processing.
#   This is also the working directory for the processing done in 
#   this script. Unzipped csv files are redirected here for processing.
#   Temporary files are stored here before being appended to the symbol
#   file csv in the appropriate directory
# filesroot/.archive
#   Directory that contains the original files after processing, sorted
#   into the following two sub-directories:
# filesroot/.archive/zip_files
# filesroot/.archive/csv_files

# Files are downloaded into yearly zip files:
# INDEX_2005.zip
# INDEX_2006.zip
# ... etc.

# Set the root of the file structure here:
# filesroot = "/home/peter/Data/EOD\ Global\ Indexes"
filesroot = "/home/peter/Data/EOD.Global.Indexes"

setwd(paste(filesroot, "/.incoming", sep=""))

start_t<-Sys.time()

# Does the directory structure exist?
if (!file.exists("../.archive")){  
  dir.create("../.archive", mode="0777")
  dir.create("../.archive/csv_files", mode="0777")
  dir.create("../.archive/zip_files", mode="0777")
}
if (!file.exists("../.archive/zip_files"))
  dir.create("../.archive/zip_files", mode="0777")

if (!file.exists("../.archive/csv_files"))
    dir.create("../.archive/csv_files", mode="0777")

# Unzip the zip files contained in filesroot/.incoming
zipfiles = list.files(pattern="*.zip")
if(length(zipfiles)>0){
  system("unzip \\*.zip")
  system("mv *.zip ../.archive/zip_files/")
} else {
  print("No zip files to process.")
}

# That creates a set of daily files, with filenames formatted as:
# INDEX_20060123.csv

# What csv files are we working with?
files = list.files(pattern="*.csv")
if(length(files) == 0){
  stop("There are no csv files to process in the .incoming directory.")
}

# Check to see if the files have already been processed
prevfiles = list.files("../.archive/csv_files")
rmfiles = files[files %in% prevfiles]
if(length(rmfiles) >0){
  file.remove(rmfiles)
  files = files[!files %in% prevfiles]
}

if (length(files) == 0)
  stop("There are no files to process or these files have been processed previously. Stopping.")
  
# Each file contains something like the following:
# Symbol,Date,Open,High,Low,Close,Volume
# ADR.IDX,23-Jan-2006,748.15,758.54,748.15,757.8,0
# ADVA.IDX,23-Jan-2006,45,559,45,549,0
# ADVN.IDX,23-Jan-2006,899,2213,899,2147,0
# ADVQ.IDX,23-Jan-2006,1440,1666,1440,1645,0
# AEX.IDX,23-Jan-2006,428.16,432.38,427.9,431.88,102237000
# AJT.IDX,23-Jan-2006,16.66,16.67,16.54,16.59,0
# ATX.IDX,23-Jan-2006,3846.97,3872.43,3804.52,3872.43,3018800
# BANK.IDX,23-Jan-2006,3088.71,3109.2,3086.72,3105.77,0
# BDI.IDX,23-Jan-2006,2417,2417,2417,2417,0

# Remove the first column and place the rest of the line in a csv file
# named for the symbol using awk
# awk -F "," 'NR!=1 {file=$1;sub($1FS,blank); print >>file ".csv"}' INDEX_20060207.csv
# Ignores the header in the original file (with NR!=1)
for (file in files){
  print(paste("Splitting ",file,sep=""))
  system(paste('awk -F "," ' , "'NR!=1 {filename=$1; sub($1FS,blank); print >> filename",'".csv"',"}'" , file, sep=" "))
  # Move the now-processed dated csv files into a 'processed' directory
  system(paste("mv ", file, " ../.archive/csv_files/", file, sep=""))
}
# That creates a .csv file for each symbol processed and moves processed
# files into the archive directory.

# Now, we want to append the resulting csv files to csv data in the
# directory for each symbol.

# What symbols do we need to process?
tmpfiles = list.files()
header = 'Date,Open,High,Low,Close,Volume'
for (file in tmpfiles){
  targetdir =  strtrim(file, (nchar(file)-4))
  fullpathdir = paste("../", targetdir, sep="")
  fullpathfile = paste("../", targetdir, "/", file, sep="")
  if (!file.exists(fullpathdir)){  # Does the directory exist?
    # No, create the directory 
    dir.create(fullpathdir, mode="0777")
    # ...and an empty file with a header within it
    system(paste("echo ", header," > ", fullpathfile, sep=""))
  }
  # Yes, directory exists
  # ... so append the local tmp file to the existing data file
  print(paste("Updating ", file, sep=""))
  system(paste("cat ", file, " >> ", fullpathfile, sep=""))
  # ... and remove the local tmp file
  file.remove(file)
}
end_t<-Sys.time()
print(c("Elapsed time: ",end_t-start_t))
print(paste("Processed ", length(files) ," days of prices for ", length(tmpfiles), " symbols.", sep=""))

# EOD provides a two-column text file, tab separated, that contains a list
# of symbols and a short description.  Its usually downloaded as
# "INDEX.txt" from somewhere on the website.  Although it is sparse, we
# need to load the list of symbols as instruments and we might as well
# keep the description (for what it's worth).

# In addition, we need to add columns for other metadata we want to 
# associate with the symbols.  Create a csv file in a spreadsheet with
# the columns from the INDEX.txt file labeled "primary_id" and "description",
# then add columns for currency (fill the column with a nonsense value
# like "EOD", for now), exchange ("EOD"), multiplier ("1"),
# source ("EODdata.com"), and any other attribute you want to associate
# with the symbol.  These suggested values are, of course, junk but
# indicate that the metadata needs to be corrected later.  With a 
# large number of symbols like this, this is just a quick way to get
# started with the data.  Once a few symbols have been selected, you 
# can correct the metadata for those contracts and use them as you 
# see fit.

# Define the nonsense currency:
# require(FinancialInstrument)
# currency(EOD)

# Then load the instruments csv file you created:
# load.instruments("~/Data/EOD.Global.Indexes/.scripts/instr.EODdata.csv")

# Now, whenever you log in you need to register the instruments.  This
# might be a line you put into .Rprofile so that it happens automatically:
# require(quantmod) # this requires a development build after revision 560 or so.
# setSymbolLookup.FI(base_dir='/home/peter/Data/EOD.Global.Indexes', split_method='common', storage_method='csv', src='csv', extension='csv', format='%d-%b-%Y')

# Now you should be able to:
# > getSymbols("FTSE.IDX")
# [1] "FTSE.IDX"
# > chart_Series(FTSE.IDX)
# > head(FTSE.IDX)




