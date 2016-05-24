`MyReportEnd` <- function#gracefully finalizes the HTML page
###gracefully finalizes the HTML page
(file = "report.html" ##<< file to append to
 ) { 
 cat("\n<hr size=1></body></html>", 
 file = file, append = TRUE) 
 }

