`readCountsData` <-
function(filename, ctPerSec)
{
    print("Please wait while I am reading your source data ...")
    #reading raw data
    Tfile <- file(filename, "r")
    if(isOpen(Tfile, "r")) #  TRUE
    {
        seek(Tfile, 0, rw="r") # reset to beginning
        lines = readLines(Tfile)
        close(Tfile)
    }

    skipPos = grep("-----", lines)[2]  #number of skip lines
    startTPos = grep("Start Time", lines)  #start time

    #get start date
    startTime = gsub("Start Time ", "", lines[startTPos])
    startTime = gsub("[[:blank:]]", "", startTime)
    startDatePos = grep("Start Date ", lines)  #startdate
    startDate = gsub("Start Date ", "", lines[startDatePos])
    startDate = gsub("[[:blank:]]", "", startDate)
    startDate = strsplit(startDate, "/")[[1]]
    if(nchar(startDate[1]) == 1){
        startDate[1] = paste("0", startDate[1], sep ="")
    }
    if(nchar(startDate[2]) == 1){
        startDate[2] = paste("0", startDate[2], sep ="")
    }
    startDate = paste(startDate[3], startDate[1], startDate[2], sep = "-")
    #end of getting startdate
    rawTimeStamp1  = paste(startDate, startTime, sep = " ")
 
    #get epochtime
    #if(is.null(ctPerSec)){
    #    ctPerSec = getEpoch(filename, unit = "sec")
    #}
    #end of getting epoch time

    startline = skipPos+1
    endline = length(lines)

    rawdata = c()
    timeline = c()
    for(i in startline: endline){
        temp0 = gsub("[[:blank:]]+", " ",lines[i])
        temp = strsplit(temp0, " ")[[1]]
        temp = temp[temp != ""]
        rawdata = c(rawdata, temp)
    }

    if(ctPerSec >1){
        timeline = rep(0:as.integer(length(rawdata)/ctPerSec), each = ctPerSec)[1:length(rawdata)]
    }else{  
        timeline = (0:as.integer(length(rawdata)-1)/ctPerSec)
    }

    rawTimeStamp = rep(rawTimeStamp1, length(rawdata))
    rst = gsub(" GMT", "", as.POSIXlt(rawTimeStamp, tz = "GMT")+ timeline)
    data.frame(TimeStamp = as.vector(rst), counts = as.numeric(as.vector(rawdata)))
}#end of readCountsData

