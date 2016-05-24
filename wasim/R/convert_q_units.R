# converts discharge data of a given file (infile) in WASIM-format into output file (outfile)
# converting discharge from [m3/s] to [mm] (convert_to="specific") or vice versa (convert_to="absolute")
# using the file gauges_area_file containing the catchment areas and a timestep of timestep_sec seconds

`convert_q_units` = function (infile, outfile, gauges_area_file, convert_to, timestep_sec)
{
#load gauge catchment sizes
gauges_area<-read.table(file= gauges_area_file,header=TRUE, sep="\t")  # read catchment sizes

infileheader=read.table(file= infile,nrows=5,sep="|",stringsAsFactors=FALSE)  # read input file header
obs_data<-read.table(file= infile,header=TRUE, na.strings =c("9999","-9999","-1.0e+04","999","-999") , skip=4)  # read input file

converted_discharge=data.frame(obs_data[,1])  #create dummy dataframe

for (i in 1:ncol(obs_data))
{
  curr_gauge=names(obs_data)[i]
  if (curr_gauge %in% c("YYYY","MM","DD","HH"))
  {
    converted_discharge[,i]=obs_data[,i]
    names(converted_discharge)[i]=curr_gauge 
    next  #skip date columns
  }
  corresponding_row=which(gauges_area$GAUGE==curr_gauge)  #get row of corresponding gauge
  if (length(corresponding_row)<1)
  {
    warning(paste("No area specified for gauge'",curr_gauge,"', runoff conversion failed."))
    converted_discharge[,i]=obs_data[,i]*NaN
    next
  }
  if (convert_to=="specific") #m3/s -> mm 
  {
    converted_discharge[,i]=obs_data[,i]*timestep_sec/(gauges_area$AREA_UPSTREAM_KM2[corresponding_row]*1e6)*1e3
  } else {         #mm -> m3/s
    converted_discharge[,i]=obs_data[,i]*(gauges_area$AREA_UPSTREAM_KM2[corresponding_row]*1e6)/1e3/timestep_sec   
  }
  names(converted_discharge)[i]=curr_gauge 
}

#summary(converted_discharge)
#summary(obs_data)



if (convert_to=="specific")
{
  newunit="[mm]"
  digits=4
 } else {
  newunit="[m3/s]"
  digits=2
}
infileheader[1,]=paste(gsub("\\[.*\\]", "", infileheader[1,]),newunit,sep="") #remove old unit, add new unit
write.table( infileheader, file= outfile, append = FALSE, quote = FALSE,row.names=FALSE,col.names=FALSE) 

na_value=-999
converted_discharge[is.na(converted_discharge)]=as.numeric(na_value)  #replace nodata values with nodata indicator
if (convert_to=="specific")
  write.table( cbind(converted_discharge[,1:4],format(round(converted_discharge[,-(1:4)],digits=digits),scientific=10)), file= outfile, append = TRUE, quote = FALSE,row.names=FALSE,col.names=FALSE,na="-999",sep="\t") else
  write.table( cbind(converted_discharge[,1:4],format(converted_discharge[,-(1:4)],digits=digits,scientific=1)), file= outfile, append = TRUE, quote = FALSE,row.names=FALSE,col.names=FALSE,na="-999",sep="\t") 

}





