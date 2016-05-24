#' Generate main window for interactive plot windows
#' 
#' \code{vlauch} generates a main window which opens each pre-generated window including statistical plot with interactivity
#' 
#' @docType methods
#' @param data data frame for default data set
#' @param name character for the name of the generated scatter plot
#' @param tag character for the common name of a series of linked plots
#' @param iframe logical; if \code{TRUE}, all plot windows are displayed in the main window as inline frames
#' @param browse logical; if \code{TRUE}, browse the main window by the default web browser through the local web server; if \code{FALSE}, generating only
#' @export
#' @examples
#' data(vsfuk2012)
#' vscat(MarriageRate, DivorceRate, vsfuk2012, "scat1", "vsfuk2012", colour=Type)
#' vhist(FertilityRate, vsfuk2012, "hist1", "vsfuk2012", fill=Type)
#' vlaunch(vsfuk2012, "main", "vsfuk2012", browse=FALSE)
#' 

vlaunch <- function(data, name, tag, iframe=FALSE, browse=TRUE){
  fn <- paste0(name, ".", tag)
  basehtmlfn <- paste0(".", tag, ".svg.html")
  htmlfn <- paste0(name, basehtmlfn)
  #fnregex <- paste0("*", basehtmlfn)
  plotfilelist <- list.files(pattern=paste0("*", basehtmlfn))
  plotfilenames <- gsub(basehtmlfn,"",plotfilelist)
  
  winlist <- paste0("var winlist=['",
                   gsub(",","','",
                        paste(plotfilenames, collapse=",")),
                   "'];\n")

  jspath <- file.path(system.file(package="vdmR"), "exec/vdmr_main.js")
  file.copy(jspath, paste0(fn, ".js"), overwrite=TRUE)
  
  csspath <- file.path(system.file(package="vdmR"), "exec/vdmr_main.css")
  file.copy(csspath, paste0(fn, ".css"), overwrite=TRUE)
  
  z <- file(paste0(fn, ".html"),"w")
  cat("<html><head><title>", file=z)
  cat(fn, file=z)
  cat("</title></head>", file=z)

  cat("<script type=\"text/javascript\"><!--\n", file=z)
  cat(winlist, file=z)
  cat(paste0("var tag='",tag,"';\n"), file=z)
  cat(paste0("var colnames= ", rjson::toJSON(colnames(data)), ";"), file=z)
  cat("\n//--></script>\n", file=z)
  
  cat("<script type=\"text/javascript\" src=\"", file=z)
  cat(paste0(fn, ".js"), file=z)
  cat("\"></script>", file=z)
  
  cat("<link rel=\"stylesheet\" type=\"text/css\" href=\"", file=z)
  cat(paste0(fn, ".css"), file=z)
  cat("\">", file=z)
  
  cat("<link rel=\"stylesheet\" type=\"text/css\" href=\"//cdn.datatables.net/1.10.0/css/jquery.dataTables.css\">", file=z)
  cat("<script type=\"text/javascript\" charset=\"utf8\" src=\"//code.jquery.com/jquery-1.10.2.min.js\"></script>", file=z)
  cat("<script type=\"text/javascript\" charset=\"utf8\" src=\"//cdn.datatables.net/1.10.0/js/jquery.dataTables.js\"></script>
", file=z)
  
  cat("<body onload=\"init();\">", file=z)

  if(iframe==FALSE){
  	cat("<div id=\"buttons\"></div><br/>", file=z)
  }

  data.ncol <- ncol(data)
  
  cat("<form id=\"form1\">Identify: <select id=\"identify\" name=\"identify\">", file=z)
  cat(paste0(rep("<option value=\"", data.ncol), 0:(data.ncol-1), "\">",
  					 colnames(data), "</option>"), file=z)
  cat("</select>", file=z)
  cat("&nbsp;&nbsp;&nbsp;Selection Box: <select id=\"selbox\" onChange=\"setSelToolVisibility(this)\"><option value=\"visible\">visible</option>", file=z)
  cat("<option value=\"hidden\">hidden</option></select></form>", file=z)

  if(iframe==TRUE){
	  cat("<div id=\"inlineplot\" style=\"width:100%; height:400\">", file=z)
  	pfnum <- length(plotfilelist)
	  for(i in 1:pfnum){
	   	cat("<iframe id=\"",plotfilenames[i],"\" scrolling=\"no\" width=\"",100/pfnum,"%\" height=\"400\" src=\"", plotfilelist[i], "\" frameborder=\"0\"></iframe>", sep="", file=z)
	  }
	  cat("</div>", file=z)
  }  
  
  cat("<table id=\"df\" class=\"display\" cellspacing=\"0\" width=\"100%\">", file=z)
  
  cat("<thead><tr>", paste(
    rep("<th>",data.ncol),
    colnames(data),
    rep("</th>",data.ncol), sep=""
  ), "</tr></thead>", sep="", file=z)
  
  cat("<tbody>", file=z)
  for(r in 1:nrow(data)){
    cat("<tr>", paste0(
      rep("<td>", data.ncol),
      as.vector(t(data[r,])),
      rep("</td>", data.ncol)
    ), "</tr>", sep="", file=z)
  }
  cat("</tbody></table>", file=z)
  
  cat("</body></html>", file=z)
  close(z)
  
  if(browse==TRUE){

    if(.Platform$GUI!='RStudio'){
      s <- Rook::Rhttpd$new()
      s$add(name="vdmR",
            app=Rook::Builder$new(
              Rook::Static$new(root=getwd(), urls="/"),
              Rook::Redirect$new(paste0("/", fn, ".html"))))
      s$start()
      s$browse(1)
    } else {
      dir <- tempfile()
      dir.create(dir)
      fcp <- file.copy(list.files(".",paste0("*.",tag,".*")), dir)
      utils::browseURL(file.path(dir,paste0(fn,".html")))
    }
  }
  
}
