### Run this example in successive copy-paste batches

## Batch 1: be careful to copy and paste only the first 3 lines 
# *without* the white-space below them.
cat("R will now open a new .xls worksheet for you!\n")
cat("Please enter path and filename, without extension:\n")
filestring<-readLines(n=1)


# R is waiting for you... enter the filename ... then proceed to next batch.

## Batch 2
newPath<-paste(filestring,'xls',sep='.')
blankbook<-XLwriteOpen(newPath)

# If you check to see whether the file exists - it's not there.
# The spreadsheet is only in R's memory. The next batch will save it.

## Batch 3
XLConnect::saveWorkbook(blankbook)
cat("Now there should be a blank file called",newPath, "- Check it out!\n")

## Now: writing into the file and resaving
# Make sure you close the file in case you opened it in Excel.
# We'll just write something silly there now:

## Batch 4
# Excel showed 1 blank sheet. But for R, there are 0 sheets until you create some.
XLConnect::createSheet(blankbook,"one") 
XLConnect::writeWorksheet(blankbook,"Something Sillee!!!",sheet='one') 
XLConnect::saveWorkbook(blankbook)

# Now it's not blank anymore - Check it out... 
# You will notice XLConnect has interpreted the string 
# as a data frame. Data transfer can only occur in the form of
# data frames (except some graphics).
# After closing the file run the last batch, which finally demonstrates 
# what XLwriteOpen itself does (open with overwrite).
# Don't forget to close the .xls file first!


## Batch 5

blankbook2<-XLwriteOpen(newPath)
XLConnect::saveWorkbook(blankbook2)

#### Now the file is blank again - Check it out!
#### All done!
