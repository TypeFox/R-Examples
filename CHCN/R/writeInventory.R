writeInventory <- function(Inventory,filename = "Inventory.inv",directory=DATA.DIRECTORY){
  
  if (!file.exists(directory)) dir.create(directory)
  fname <- file.path(directory,filename,fsep =.Platform$file.sep)
  columnOrder <- c(6,3,4,5,1,2,7,8,9)
  write.table(Inventory[,columnOrder],fname )
}