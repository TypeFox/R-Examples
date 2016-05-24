downloadRtopExampleData = function(folder = system.file("extdata",package="rtop")) {
  wd = getwd()
  setwd(folder)
  download.file("http://www.hydro.tuwien.ac.at/fileadmin/mediapool-hydro/Downloads/rtopData.zip", 
         "rtopData.zip")
  unzip("rtopData.zip")
  setwd(wd)
}