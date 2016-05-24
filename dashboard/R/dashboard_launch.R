#' dashboard_launch generates the dashboard and launchs the local server
#'
#' \code{dashboard_launch} writes all files for the web page and launchs the local server
#' 
#' @docType methods
#' @param browse boolean for launching the local server
#' @export
#' @examples
#' dashboard_open(data=iris) # other options: pathoutput=getwd() ...
#' dcpiechart(x=names(iris)[5])
#' dcbarchart(x=names(iris)[1] , gap=75)
#' dcpiechart(x=names(iris)[2])
#' dctable(index=names(iris)[5])
#' dashboard_launch(browse = FALSE) # Just generates files. Server is not launched
#'  

dashboard_launch <-function(browse = TRUE){

# internal function : copy all library required in the ouput folder
.copyfilename <- function(fn){
	fnpath <- file.path(system.file(package="dashboard"), paste("exec/",fn, sep=""))
	file.copy(fnpath, paste(dashboard.env$path_folder_output,"/lib/", fn, sep=""), overwrite=TRUE)	
	Sys.chmod(paste(dashboard.env$path_folder_output, "/lib/", fn, sep=""), "777")
}

# copy all library required in the ouput folder	
.copyfilename("bootstrap.min.css")
.copyfilename("bootstrap.min.js")
.copyfilename("dc.css")
.copyfilename("dc.min.js")
.copyfilename("jquery-1.11.1.min.js")
.copyfilename("crossfilter.min.js")
.copyfilename("d3.min.js")
.copyfilename("jquery.dynatable.css")
.copyfilename("jquery.dynatable.js")

	# Write all in output_file

# for html

begin <- paste("<!DOCTYPE html><html lang='en'>
	<head>
	<meta charset='utf-8'>
	<title>",dashboard.env$title,"</title>
	<script src='lib/d3.min.js' type='text/javascript'></script>
	<script src='lib/crossfilter.min.js' type='text/javascript'>
	</script><script src='lib/dc.min.js' type='text/javascript'>
	</script><script src='lib/jquery-1.11.1.min.js' type='text/javascript'>
	</script><script src='lib/bootstrap.min.js' type='text/javascript'>
	</script><script src='lib/jquery.dynatable.js' type='text/javascript'>
	</script><link href='lib/jquery.dynatable.css' rel='stylesheet' type='text/css'>
	<link href='lib/bootstrap.min.css' rel='stylesheet' type='text/css'>
	<link href='lib/dc.css' rel='stylesheet' type='text/css'>
	<style>h4 span {font-size:14px; font-weight:normal;}</style>
	</head>
	<body>
	<div class='container' style='font: 12px sans-serif;'>
		<div class='row'><h2><div  style='float: left;'>",dashboard.env$title,"<br>
		<span  class='dc-data-count'><span class='filter-count'></span> &nbsp;<span class='filter-count-text'>selected out of</span>&nbsp;<span class='total-count'></span>&nbsp;records <a class='dc-data-count-reset-all' href='javascript:dc.filterAll(); dc.renderAll();'>Reset All</a></span>
		</div></h2></div>
		<div class='row'><div class='span12'><h2> &nbsp;</h2></div></div>
		<div class='row'>",sep="")

end <- paste("</div><div class='row'><div class='span12'></div></div>
	</div><script src='",dashboard.env$filename,".js' type='text/javascript'></script>
	</body></html>",sep="")

Sys.chmod(paste(dashboard.env$path_folder_output, "/lib", sep=""), "777")

Sys.chmod(paste(dashboard.env$path_folder_output, "/",dashboard.env$filename, ".csv", sep=""), "777")

# write the html file
fhtml <- file(paste(dashboard.env$path_folder_output, "/",dashboard.env$filename, ".html", sep=""),"w")
  cat(begin, file=fhtml)
  cat(dashboard.env$temp, file=fhtml)
  cat(end, file=fhtml)
  close(fhtml)

Sys.chmod(paste(dashboard.env$path_folder_output, "/",dashboard.env$filename, ".html", sep=""), "777")

# for javascript

# write in the js file : 
fjs <- file(paste(dashboard.env$path_folder_output, "/",dashboard.env$filename, ".js", sep=""),"w")
  cat(dashboard.env$jstemphead , file=fjs)
  cat(dashboard.env$jstempcore, file=fjs)
  cat(dashboard.env$jstempendt1, file=fjs)
  cat(dashboard.env$jstempend, file=fjs)
  cat(dashboard.env$jstempendt3, file=fjs)
  cat(" dc.renderAll();}", file=fjs)
  close(fjs)

Sys.chmod(paste(dashboard.env$path_folder_output, "/",dashboard.env$filename, ".js", sep=""), "777")

Sys.chmod(dashboard.env$path_folder_output, "777")


	# TODO: Export all in exportzip (all js lib would be in the html file)

	# start the server and open the dashboard web page
  if(browse==TRUE){
    dashboard.env$s <- Rook::Rhttpd$new()
    dashboard.env$s$add(name="dashboard",
          app=Rook::Builder$new(
            Rook::Static$new(
              root=dashboard.env$pathoutput,
              urls="/"),
            Rook::Redirect$new(paste("/", dashboard.env$folder_output , "/", dashboard.env$filename, ".html", sep=""))))
    dashboard.env$s$start()
    dashboard.env$s$browse(1)
   }
  
  cat("Files have been generated in this folder: ",dashboard.env$path_folder_output)	
	
}



 