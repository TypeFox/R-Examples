### Name: read.observations
### Title: Reads runoff-observation in WASIM-format,trim to extent of simulated data
### Aliases: read.observations
### Keywords: ~load ~read ~import ~WASIM data ~observation ~observed ~measured discharge

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

## The function is currently defined as
`read.observations` <-
function (filelist,       #list of input files to be considered
  path="wasim/input",     #path where to look for input files
  gauge_names=NULL,       #names of gauges that should be read (NULL:read all)
  date_modelled=NULL,          #timespan that will be extracted
  na.strings=c("999", "999.00","9999", "9999.00", "9999.000000","-9999","-9999.0","-9999.00"), #strings indicating missing data
  sep="" #field separator in input file
)
{

read.files.nr=length(filelist)    #get number of files to be read
if ((nchar(path)>0) && (path[length(path)]!="/") && (path[length(path)]!="\\"))    #remove trailing (back)slash
  path = paste(path, "/",sep="")

if(!is.null(date_modelled)){
     obs_data_all=data.frame(datenum=date_modelled)
     all.y=FALSE
} else {
     obs_data_all=NULL
     all.y=TRUE
}
for (i in 1:read.files.nr ){
 filename <- paste(sep="", path, filelist[i])
 cat("i: ", i, " ", filename, "\n")
 curr_gauge_names<-read.table(file=filename,header=FALSE, skip=4,nrows=1,sep=sep,stringsAsFactors=FALSE)  # retrieve  column names
 
  if (is.null(gauge_names))
    relevant_columns=5:ncol(curr_gauge_names) else  #read all
    relevant_columns=which(curr_gauge_names[1,] %in% gauge_names)  #look for desired gauge names
    
 if (length(relevant_columns)>0)
 {
      colClasses=rep("NULL",ncol(curr_gauge_names))      #skip all columns...
      colClasses[c(1:4,relevant_columns)]="numeric"           #...except these

      obs_data = read.table(filename[i],row.names=NULL, na.strings=na.strings, skip=4,colClasses=colClasses,header=TRUE,sep=sep)
      names(obs_data)[1:4] = c("Year", "Month", "Day", "Hour")
      obs_data$datenum=as.POSIXct(assemble.date(obs_data[,1:4]),tz="GMT")

      
      if(is.null(obs_data_all)){
          obs_data_all=obs_data[,-(1:4)]
      } else {
          obs_data_all=merge(obs_data_all,obs_data[,-(1:4)],by="datenum",all.x=TRUE, all.y=all.y)
      }
   
      gauge_names=gauge_names[-which(gauge_names %in% curr_gauge_names)]    #remove loaded gauge names from to-do-list
      if (is.null(gauge_names))
        break   #all necessary gauges have been read
        
 }

}

return(obs_data_all)

}

