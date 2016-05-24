`BasicHTML` <- structure(function#creates a basic HTML page displaying plots and annota
###   Creates a basic HTML page displaying plots and annotations that can easily be navigated. The plots can be created either 'on the fly' by passing the appropriate commands or beforehand in which case just the filenames need to be passed.                   
( cmds=NULL,  ##<< list of commands that generates the plots. If missing, the graphfiles are assumed to exist already.
  HTMLobjects,  ##<< list of graph filenames, either to be created by the list of commands or to be copied to the Figures subdirectory and/or dataframes to be displayed in sortable tables.
  Captions,  ##<< vector of captions; these go directly below the graphs
  MenuLabels,  ##<< vector of labels for the main page.
  Comments=NULL,  ##<< Text/comments to be written between the graphs
  file="tmp.html", ##<< file name of main page; '.html' extension will be added. The '\_main' and '\_menu' pages use this base as well.
  title="", ##<< title to be written in the navigation/menu page
  width=480, ##<< width for all graphfiles
  height=480, ##<< height for all graphfiles
  FRAMES=FALSE, ##<< is this an HTML page with frames ?
  JSCPATH = "jsc", ##<< path that should contain the jsc components. If non existing, user will be prompted for installation.
  LaunchPage = FALSE, ##<< launch the page ?
  APPEND = FALSE,   ##<< append to existing HTML page ?
  href = NULL,  ##<< links to other HTML pages
  verbose=0 ##<< level of verbosity
){
	#BasicHTML(cmds = c("plot(rnorm(100));","plot(1:10);"), HTMLobjects = list("Fig1.png", "Fig2.png"), Captions=c("Gaussian noise","seq 1:10"))
	#x <- cbind.data.frame(x1 = round(rnorm(10),3), x2 = round(runif(10),3));
	#attr(x, "HEADER") <- "some random numbers";
	#BasicHTML(HTMLobjects = list("Fig1.png", x, "Fig2.png"), Captions=c("Gaussian noise","Gaussian and uniform random numbers", "seq 1:10"), file = paste(Sys.getenv("HOME"), "/public_html/tmp/tmp.html",sep=""), JSCPATH = "../jsc")
	
  HTMLInsertGraph <- function (GraphFileName = "", Caption = "", GraphBorder = 1, 
    Align = "center", WidthHTML = 500, HeightHTML = NULL, file = get(".HTML.file"), 
    append = TRUE, href = NULL, ...) 
   {
     cat("\n", file = file, append = append, ...)
    
    cat(paste("<p align=", Align, ">",if (!is.null(href)) 
            paste("<a href=\"", href, "\"> ", sep = "")
        else "","<img src='", GraphFileName, 
        "' border=", GraphBorder, if (!is.null(WidthHTML)) 
            paste(" width=", WidthHTML, sep = "")
        else "", if (!is.null(HeightHTML)) 
            paste(" height=", HeightHTML, sep = "")
        else "", ">", if (!is.null(href)) "</a>" else "", sep = "", collapse = ""), file = file, 
        append = TRUE, sep = "")
    if (Caption != "") 
        cat(paste("<br><i class=caption>", Caption, "</i>"), 
            file = file, append = TRUE, sep = "")
    invisible(return(TRUE))
   }

  #require(Cairo);
  graphIndex <- sapply(HTMLobjects,is.character);
  graphfiles <- unlist(HTMLobjects[graphIndex]);
  header <- NULL;
  if (!all(graphIndex)) {
  	header <- HTMLsortedTable(x=NULL, JSCPATH = JSCPATH);
  	if (!file.exists(paste(JSCPATH,sep="")))  InstallJSC(paste(JSCPATH,sep=""));
  }
  
  if (is.null(cmds)){#no cmds passed, so we assume that the plots have already been created !
  	DoNotPlot <- TRUE;
  } else  {
  	stopifnot(length(cmds)==length(graphfiles));
    DoNotPlot <- FALSE;
    #CairoWorks();
  }
  
  if (verbose >1) browser();
  
  if(! FRAMES & !APPEND ) MyReportBegin(file=file, title = title, header = header);
  if (APPEND & FRAMES ){
  	tmp <- scan(paste(substring(file,1,nchar(file)-10), "_menu.html",sep=""), what = "");
  	ExistingLabels <- grep("#Num", tmp);
  	MenuNumber <- length(ExistingLabels) + 1;
  } else {
    MenuNumber <- 1;}
  
  for (i in seq(along= HTMLobjects)){
  	if (FRAMES){
  	  if (!is.null(Comments)) if (nchar(Comments[i])>0) HTML(Comments[i]);
  	  MenuLabel <- ""; 
  	  if (!missing(MenuLabels)) if (!is.null(MenuLabels)) if (length(MenuLabels) == length(HTMLobjects)) {MenuLabel <- MenuLabels[i]};
  	  HTML(paste("<a name=Num", MenuNumber,">&nbsp;</a><p><xmp class=command>", MenuLabel,"</xmp>",sep=""),file= file)
      MenuNumber = MenuNumber + 1;
    }#end of if (FRAMES)  
    #if (is.character(HTMLobjects[[i]]) ){
    if (graphIndex[i]){
      graphfile <- HTMLobjects[[i]];
  	  if (!DoNotPlot){
        png(graphfile,width=width,height=height);
          if (verbose) cat("creating graph", graphfile,"\n")
	        eval(parse(text=cmds[i]));
	      if (!is.null(dev.list())) dev.off()
	  }
	  HTMLInsertGraph(graphfile,Caption=Captions[i], Align="left",file=file,WidthHTML=width, HeightHTML=height, href=href[i]);
	} else if (is.data.frame(HTMLobjects[[i]]) | is.matrix(HTMLobjects[[i]])) {
	  HEADER <- "";
	  try(HEADER <- attr(HTMLobjects[[i]],"HEADER"));
	  htable <- HTMLsortedTable(HTMLobjects[[i]], HEADER = HEADER, file = NULL, JSCPATH = JSCPATH);
	  HTML('<p align= center >', file = file);HTML(htable, file = file);HTML('<BR>', file = file);
	  HTML(Captions[i], file = file);HTML('<BR>', file = file);
	  #HTML(x, file = file)
	}
  }


 #for now I am not adding a clean end of html footer in order to keep the pages open for appending operations
 #most browsers seem to be able to handle this just fine...
 # if(! FRAMES) MyReportEnd(file=file);
  
  ##seealso<< \link{FramedHTML}
  
  if (LaunchPage) system(paste("open ", file))
### no return value
}, ex=function(){
  if (interactive()){  
  owd=setwd(tempdir())
  BasicHTML(cmds = list("plot(rnorm(100));","plot(1:10);"), 
            HTMLobjects = list("Fig1.png", "Fig2.png"), 
            Captions=c("Gaussian noise","seq 1:10"),  
            MenuLabels = c("Marvel at the graph below","scatterplots are nice"), 
            title="Test Page",width=480, height=480, verbose=1, JSCPATH = NULL)
                      
    #example with plots and graphfiles having been generated beforehand:
    png("Fig1.png");
      plot(rnorm(100));
    dev.off()
    png("Fig2.png");
      plot(1:10);
    dev.off();
    
BasicHTML( HTMLobjects = list("Fig1.png", "Fig2.png"), 
  Captions=c("Gaussian noise","seq 1:10"), 
   MenuLabels = c("Marvel at the graph below","scatterplots are nice"), 
  title="Test Page",
  width=480, height=480, verbose=1, JSCPATH = NULL);
    
    #example with absolute paths for graphfiles :
    Fig1 <- paste(tempdir(),"/Fig1.png",sep="")
    png(Fig1);
      plot(rnorm(100));
    dev.off()
    Fig2 <- paste(tempdir(),"/Fig2.png",sep="")
    png(Fig2);
      plot(1:10);
    dev.off();
    
BasicHTML( HTMLobjects = list(Fig1, Fig2), 
    Captions=c("Gaussian noise","seq 1:10"),  
    MenuLabels = c("Marvel at the graph below","scatterplots are nice"), title="Test Page",
    width=480, height=480, verbose=1, JSCPATH = NULL);
    #cleanup:
    #system(paste("rm ", Fig1));system(paste("rm ", Fig2))

  #example with sorted table:
  x <- cbind.data.frame(x1 = round(rnorm(10),3), x2 = round(runif(10),3));
  attr(x, "HEADER") <- "some random numbers";
  BasicHTML(HTMLobjects = list("Fig1.png", x, "Fig2.png"), 
            Captions=c("Gaussian noise","Gaussian and uniform random numbers", "seq 1:10"), 
            file = paste(Sys.getenv("HOME"), "/public_html/tmp/tmp.html",sep=""), 
            JSCPATH = "../jsc");
  setwd(owd)
}
})
