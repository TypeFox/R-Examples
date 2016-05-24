testSWAT2005=function(){
runSWAT2005=function (hist_wx = NULL, elev = 100, rch = 3) {
     Sys.setenv(GFORTRAN_STDIN_UNIT = -1)
     Sys.setenv(GFORTRAN_STDOUT_UNIT = -1)
     Sys.setenv(GFORTRAN_STDERR_UNIT = -1)
     if (length(hist_wx) > 2) {
         tmp_head = paste("Tmp\nLati Not Used\nLong Not Used\nElev        ", 
             sprintf("%5.0f\n", elev), sep = "")
         pcp_head = paste("Pcp\nLati Not Used\nLong Not Used\nElev        ", 
             sprintf("%5.0f\n", elev), sep = "")
         cat(tmp_head, sprintf("%s%005.1f%005.1f\n", format(hist_wx$DATE, 
             "%Y%j"), hist_wx$TMX, hist_wx$TMN), file = "tmp.tmp", 
             sep = "")
         cat(pcp_head, sprintf("%s%005.1f\n", format(hist_wx$DATE, 
             "%Y%j"), hist_wx$PRECIP), file = "pcp.pcp", sep = "")
         print("built new pcp.pcp and tmp.tmp files, make sure they are correct in file.cio")
     }
     libarch = if (nzchar(version$arch)) 
         paste("libs", version$arch, sep = "/")
     else "libs"
     swatbin <- "rswat2005.exe"
     system(shQuote(paste(path.package("SWATmodel"), libarch, swatbin, sep = "/")))
     start_year = read.fortran(textConnection(readLines("file.cio")[9]), 
         "f20")
     temp = readLines(file("output.rch"))
     rchcolname = sub(" ", "", (substr(temp[9], 50, 61)))
     flow = data.frame(as.numeric(as.character(substr(temp[10:length(temp)], 
         50, 61))))
     colnames(flow) = rchcolname
     reach = data.frame(as.numeric(as.character(substr(temp[10:length(temp)], 
         8, 10))))
     rchcolname = sub(" ", "", (substr(temp[9], 8, 10)))
     colnames(reach) = rchcolname
     outdata = cbind(reach, flow)
     temp2 = subset(outdata, outdata$RCH == rch)
     temp2$mdate = as.Date(row(temp2)[, 1], origin = paste(start_year - 
         1, "-12-31", sep = ""))
     return(temp2)
 }
dirname="testswat"
swat_general <- NULL
rm(swat_general)
load(paste(path.package("SWATmodel"),"data/swat_general.rda", sep = "/"))
dir.create(dirname)
setwd(dirname)
for (file in names(swat_general)) {
   print(file)
   cat(unlist(swat_general[file]), file = file, sep = "\n")
}

test=runSWAT2005()
if(sum(test$FLOW_OUTcms) < 200){
  print("Your computer is having SWAT array issues, please contact drf28@cornell.edu as we would like to figure why this is happening.\n")
  } else {
  print("Your computer CPU test out OK for SWAT2005.\n")
  }
  
file.remove(list.files("../testswat/"))
setwd("../")
file.remove("testswat/")

}
