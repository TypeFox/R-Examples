downloadMaster()
Stations <- writeMonthlyStations()
scrapeToCsv(Stations)
MISSING <- getMissingScrape()
EMPTY <- getEmptyCsv()
if (is.null(EMPTY) & is.null(MISSING)){
  data <- createDataset()
  ###  save all the data
  writeData(data)
  inv  <- createInventory()
  write.table(inv,"masterInventory.inv")
  # write a ghcn style inventory
  writeInventory(inv)
  # select Tave data
  Mean <- formatGhcn(data, dataColumn = 7)
  writeGhcn(Mean)
 
}  
 