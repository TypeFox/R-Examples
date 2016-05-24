`HTMLsortedTable` <-function#create sortable table 
### create sortable table using JavaScript components in \code{JSCPATH} directory
(x, ##<< data frame or matrix with column names 
 TITLE="", ##<< title for the HTML page
 HEADER = "", ##<< header to display for the sorted table
 file="tmp.html", ##<< file name of main page; '.html' extension will be added. The '\_main' and '\_menu' pages use this base as well.
 JSCPATH = "jsc", ##<< path that should contain the jsc components. If non existing, user will be prompted for installation. 
 path = paste(Sys.getenv("HOME"), "/public_html/",sep=""),##<<  directory to create the file in
 debug=0 ##<<level of verbosity
){
#left to do: specify the paths !! Right now, this works only in the main public_html directory that contains the jsc/ folder!
#example: x <- cbind.data.frame(x1 = round(rnorm(10),3), x2 = round(runif(10),3));
#         HTMLsortedTable(x, TITLE= "some random numbers", file = NULL)
#         HTMLsortedTable(x, TITLE= "some random numbers",path = paste(Sys.getenv("HOME"), "/public_html/metroPCS/",sep=""), file = "tmp.html", JSCPATH = "../jsc")
header0 <- paste('<html xmlns=\"http://www.w3.org/1999/xhtml\" \n  xml:lang=\"en\"> \n');
header <- paste(' <head> \n    <title>',TITLE,'</title>
  <link rel=\"stylesheet\" type=\"text/css\" href=\"',JSCPATH,'/jsComponents.css\" />
  <script type=\"text/javascript\" src=\"',JSCPATH,'/extra/jsRegionChart.js\"> </script>
  <script type=\"text/javascript\" src=\"',JSCPATH,'/extra/jsBarChart.js\"> </script>

  <script type=\"text/javascript\" src=\"',JSCPATH,'/jsComponents.js\"> </script>
  <script type=\"text/javascript\" src=\"',JSCPATH,'/extra/jsDebugger.js\"></script>
  <script type=\"text/javascript\" src=\"',JSCPATH,'/extra/jsDebuggerTools.js\"></script>
  <script type=\"text/javascript\" src=\"',JSCPATH,'/extra/jsHelloWorld.js\"></script>
  <link rel=\"stylesheet\" type=\"text/css\" href=\"',JSCPATH,'/extra/jsDebugger.css\" />
  <link rel=\"stylesheet\" type=\"text/css\" href=\"',JSCPATH,'/extra/jsHelloWorld.css\" />
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\" />
<meta http-equiv=\"Content-Style-Type\" content=\"text/css\" />
</head>',sep="");

if (is.null(x)) return(header);

#<body>
if (!is.null(HEADER))
  if (nchar(HEADER) > 0) HEADER <- paste('<h2>',HEADER,'</h2> \n');
tableheader <- paste(HEADER,' <div class=\"JSTableStripe\"> \n <div class=\"JSTableSort\"> \n <table>\n<thead>\n <tr>',sep="");

tablefooter <- paste('\n</table>\n</div></div>');
footer <- paste('</body>
</html>');

  #tableheader <- paste('<thead>\n   <tr>');
  x.col.class <- sapply(x,"class");
  ColClasses <- rep('<th class = "SortString">',ncol(x));
  ColClasses[x.col.class == "numeric"] <- '<th class = "SortNumber">';  
   	
  for (i in 1:ncol(x)) tableheader <- c(tableheader,paste('\n\t\t', ColClasses[i], colnames(x)[i], '</th>', sep = "", collapse = ""))

  tableheader <- c(tableheader, paste('</tr>\n    </thead> \n'));
  
  tablebody <- paste('<tbody> \n')
  for (j in 1:nrow(x)){
  	tablebody <- c(tablebody, paste('<tr>'));
    for (i in 1:ncol(x)) tablebody <- c(tablebody,paste('<td>', x[j,i], '</td>', sep = "", collapse = ""));
    tablebody <- c(tablebody, paste('</tr> \n'));
  }
  if (debug) cat(header, tableheader, tablebody , footer);
  if (!is.null(file)) {
  	if (!file.exists(paste(path, JSCPATH,sep="")))  InstallJSC(paste(path, JSCPATH,sep=""));
    cat(header0, header, tableheader, tablebody , footer, file = paste(path,file,sep=""));
  } else {
  	return( c(tableheader, tablebody, tablefooter));
  	#return(list(header=header, table = c(tableheader, tablebody)));
  }
  
 }

