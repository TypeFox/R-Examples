warcolak.convert <-
function(w){
# warcolak.convert() - convert warcolak data file to format required for a 
#                      dataframe for dmm() or mdf()
  df <- w
  colnames(df) <- c("Id","SId","DId","Sex","Trait1","Trait2",colnames(w)[7:13])
  df[,"Id"] <- w[,"ID"]
  df[,"SId"] <- w[,"Sire"]
  df[,"DId"] <- w[,"Dam"]
  df[,"Sex"] <- w[,"sex"]
  df[,"Trait1"] <- w[,"trait1"]
  df[,"Trait2"] <- w[,"trait2"]
  return(df)
}
