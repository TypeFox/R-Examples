`FramedHTML` <- structure(function# creates a framed HTML page displaying plots and annotations
###  Creates a framed HTML page displaying plots and annotations that can easily be navigated. The plots can be created either 'on the fly' by passing the appropriate commands or beforehand in which case just the filenames need to be passed.
### The user has a great deal of flexibility in choosing appropriate directory structures.
(  cmds=NULL, ##<< list of commands that generates the plots. If missing, the graphfiles are assumed to exist already.
   basepath = c("./", paste(Sys.getenv("HOME"), "/public_html/",sep=""))[1],  ##<< base path of \samp{public\_html} directory
   path="tmp",  ##<< subdirectory of \code{basepath}; will be created if non existing
   Graphpath = "Figures/", DiagnosticsPath = "Diagnostics",  ##<< subdirectory of \samp{basepath/path/} containing the graphfiles; will be created if non existing
   file="tmp",   ##<< file name of main page; '.html' extension will be added. The '\_main' and '\_menu' pages use this base as well.
   HTMLobjects,  ##<< list of graph filenames, either to be created by the list of commands or to be copied to the Figures subdirectory and/or dataframes to be displayed in sortable tables.
   Captions,  ##<< vector of captions; these go directly below the graphs
   MenuLabels1,  ##<< vector of labels for the menu navigation page. It helps to keep these succinct and short !.
   MenuLabels2,  ##<< vector of labels for the main page; these go on top of the individual graphs, so they are complementary to the captions.
   href = NULL,##<< links to other HTML pages
   Comments = NULL,  ##<< Text/comments to be written between the graphs
   title="", ##<< title to be written in the navigation/menu page
   width=480, ##<< width for all graphfiles
   height=480, ##<< height for all graphfiles
   FRAMES=FALSE, ##<< is this an HTML page with frames ?
   JSCPATH = "jsc", ##<< path that should contain the jsc components. If non existing, user will be prompted for installation.
   REFRESH = "",##<< Meta refresh is a method of instructing a web browser to automatically refresh the current web page after a given time interval
   img.logo.path =paste(Sys.getenv("HOME"), "/public_html/",sep=""),  ##<< path to search for the logo pic in the frame
   img.logo = NULL, ##<< filename of logo to display
   img.href= 'http://www.sensenetworks.com',  ##<< link of logo to point to.
   APPEND = FALSE,   ##<< append to existing HTML page ?
   verbose=1 ##<< level of verbosity
  ){
  #FramedHTML(cmds = c("plot(rnorm(100));","plot(1:10);"), HTMLobjects =list("Fig1.png", "Fig2.png"), Captions=c("Gaussian noise","seq 1:10"), MenuLabels1 = c("Label1","Label2"), MenuLabels2 = c("Marvel at the graph below","scatterplots are nice"), title="Test Page",verbose=1);
  #x <- cbind.data.frame(x1 = round(rnorm(10),3), x2 = round(runif(10),3));
  #FramedHTML(HTMLobjects =list(x,"Fig1.png", "Fig2.png"), Captions=c("jadejade","Gaussian noise","seq 1:10"), MenuLabels1 = c("sorted table","Label1","Label2"), MenuLabels2 = c("the magic of javascript","Marvel at the graph below","scatterplots are nice"), title="Test Page",verbose=1)

  ##note<<  There is not much eror checking. In particular, the lengths of the arguments\code{cmds, graphfiles, Captions, MenuLabels1, MenuLabels2} need to be all the same !

  ##seealso<< \link{BasicHTML} 
  
  #ArgNames <- names(formals(FramedHTML));
  #Args <- match.arg(ArgNames );
  #Args <- list();
  #return(Args)
   
  WD <- getwd();
  on.exit(setwd(WD))
  #The paths can be confusing as
  # (i)  the HTML file only wants relative paths, but Cairo wants absolute.
  # (ii) we don't know if the user passes the graphfiles with full, partial or no paths !
  #Solution: change current dir to the base path and assume that the graphfiles are just file names, no path info
  basepath <- makePathName(basepath, TRUE);
  setwd(basepath);
  path <- makePathName(path, TRUE);
  NoDirs <- sum("/" == unlist(strsplit(path,"")));#very fragile method of determining the number of subdirectories off the basepath !!
  setwd(path);
  Graphpath <- makePathName(Graphpath, TRUE);
   #If string, assume graphfile:
  graphIndex <- sapply(HTMLobjects,is.character);
  graphfiles <- unlist(HTMLobjects[graphIndex]);
 
  if (verbose > 1) browser();
  
  if (!is.null(graphfiles)){
    for (g in graphfiles) {
    	graphdir <- dirname(g);
    	gg <- basename(g);
    	if (graphdir == "." ) {ggg <- paste(WD, gg,sep="/");#no graphfile paths
    	} else if (substring(graphdir,1,1) != "/") {ggg <- paste(WD, g,sep="/");#relative graphfile paths 
    	} else {ggg <- g}#absolute paths
    	if (verbose) print(ggg);
      #if (file.exists(ggg) & !file.exists(paste(Graphpath, g,sep="/"))) {
      if (file.exists(ggg) ) {
      	if (verbose) cat("copying ", ggg, " to ", paste(basepath, path, Graphpath, gg, sep=""), "\n");
      	#system(paste("cp ", paste(WD, g,sep="/"), " ", Graphpath, g,sep="") );
      	system(paste("cp ", ggg, " ", paste(basepath, path, Graphpath, gg, sep="") ,sep="") );
      } else {cat(ggg, " does not exist, unable to copy !\n");}
    }
    graphfiles <- paste(Graphpath, basename(graphfiles),sep="");
    HTMLobjects[graphIndex] <- graphfiles;
    if (verbose) print(graphfiles)
  }
  DiagnosticsPath <- makePathName(DiagnosticsPath, TRUE);
  #FullGraphpath <- makePathName(paste(path, Graphpath, sep =""), TRUE);
  #DiagnosticsPath <- makePathName(paste(path,DiagnosticsPath,sep=""), TRUE);
   JSCPATH <- paste(paste(rep("../", NoDirs),collapse=""), "jsc",sep="");
   if (substring(file,1,1) != "/") {outdir = "."} else {outdir = ""}
  
  targetBig <- myHTMLInitFile(outdir = outdir, file, HTMLframe =TRUE, NavTitle = title, 
                              Title = "", JSCPATH= JSCPATH, useLaTeX = FALSE, REFRESH = REFRESH, 
                              APPEND = APPEND, img.logo.path =img.logo.path, img.logo = img.logo, img.href=img.href)
  
  target <- targetBig["target"]
  target.menu <- targetBig["targetmenu"]
  target.main <- targetBig["targetmain"]
  file = target.main;
  HTMLCSS(file = file, CSSfile = "R2HTML");
  #HTML.title(as.title(title), HR=2, file=file);
  if (!APPEND) HTML.title(as.title(title), HR=2, file=target.menu);
  
  if (APPEND){
  	tmp <- scan(target.menu, what = "");
  	ExistingLabels <- grep("#Num", tmp);
  	MenuNumber <- length(ExistingLabels) + 1;
  } else {
    MenuNumber <- 1;
  }
    
    if (verbose >1) {
    	print(HTMLobjects);
    	print(targetBig);
    	print(getwd())
    }
  BasicHTML(cmds, HTMLobjects, Captions, MenuLabels2, href=href, Comments, file, title, width=width, height, FRAMES=TRUE, JSCPATH= JSCPATH, APPEND = APPEND, verbose=verbose);

  for (i in seq(along= HTMLobjects)){
  	MenuLabel <- ""; 
  	  if (!missing(MenuLabels1)) if (!is.null(MenuLabels1)) if (length(MenuLabels1) == length(HTMLobjects)) {MenuLabel <- MenuLabels1[i]};
  	tmp <- paste("href='",target.main,"#Num", MenuNumber,"'",sep="");
  	tt <- paste("<a class=command ", tmp," target=main> ", MenuLabel," </a>",sep="");
    HTMLli(tt,file= target.menu);
    cat('\n', file = target.menu, append = TRUE)
    if (verbose >1) print(tt)
    MenuNumber = MenuNumber + 1;
  }
  

  graphics.off();
  
  #for now I am not adding a clean end of html footer in order to keep the pages open for appending operations
 #most browsers seem to be able to handle this just fine...
 # MyReportEnd(file=target.main);
  if (!is.null(WD)) setwd(WD);
  #return(targetBig)
### no return values
}, ex=function(){
if (interactive()){  
  #example with plots and graphfiles being generated on the fly:
  owd=setwd(tempdir())
  system("mkdir Figures")
  
FramedHTML(cmds = list("plot(rnorm(100));","plot(1:10);"), 
           HTMLobjects =list("Fig1.png", "Fig2.png"), 
           Captions=c("Gaussian noise","seq 1:10"),
           MenuLabels1 = c("Label1","Label2"), 
           MenuLabels2 = c("Marvel at the graph below","scatterplots are nice"), 
           Comments  = c("100 random numbers","Simple plot"), title="Test Page",
           width=480, height=480, verbose=1)
    
    
    #example with plots and graphfiles having been generated beforehand:
    png("Fig1.png");
      plot(rnorm(100));
    dev.off()
    png("Fig2.png");
      plot(1:10);
    dev.off();
    
FramedHTML( HTMLobjects = list("Fig1.png", "Fig2.png"), 
   Captions=c("Gaussian noise","seq 1:10"), 
  MenuLabels1 = c("Label1","Label2"), 
   MenuLabels2 = c("Marvel at the graph below","scatterplots are nice"), 
   Comments  = c("100 random numbers","Simple plot"), title="Test Page",
  width=480, height=480, verbose=1);
    
    #example with absolute paths for graphfiles :
    Fig1 <- paste(tempdir(),"/Fig1.png",sep="")
    png(Fig1);
      plot(rnorm(100));
    dev.off()
    Fig2 <- paste(tempdir(),"/Fig2.png",sep="")
    png(Fig2);
      plot(1:10);
    dev.off();
    
 FramedHTML( HTMLobjects = list(Fig1, Fig2), Captions=c("Gaussian noise","seq 1:10"), 
    MenuLabels1 = c("Label1","Label2"), 
    MenuLabels2 = c("Marvel at the graph below","scatterplots are nice"), 
    Comments  = c("100 random numbers","Simple plot"), 
    title="Test Page",width=480, height=480, verbose=1);
    #cleanup:
    #system(paste("rm ", Fig1));system(paste("rm ", Fig2))
  
  #example with sorted table:
  x <- cbind.data.frame(x1 = round(rnorm(10),3), x2 = round(runif(10),3));
  attr(x, "HEADER") <- "some random numbers";
  FramedHTML(HTMLobjects = list("Fig1.png", x, "Fig2.png"), 
    MenuLabels1 = c("Label1","Label2","Label3"), 
    MenuLabels2 = c("Marvel at the graph below","JavaScript rocks","scatterplots are nice"),
    Captions=c("Gaussian noise","Gaussian and uniform random numbers", "seq 1:10"),Comments = NULL,
    path = "tmp", file = "index");

  #example with sorted tables only, no figures:
  x <- cbind.data.frame(x1 = round(rnorm(10),3), x2 = round(runif(10),3));
  attr(x, "HEADER") <- "some random numbers";
  y <- cbind.data.frame(y1 = rbinom(10,50,0.3), y2 = rbinom(10,100,0.15));
  attr(y, "HEADER") <- "rbinom";
  FramedHTML(HTMLobjects = list( x, y), 
           MenuLabels1 = c("x","y"), 
           MenuLabels2 = c("JavaScript rocks","Secret numbers"),
           Captions=c("Gaussian and uniform random numbers", "Binomial draws"),Comments = NULL,
           path = "tmp", file = "index");

  setwd(owd)
}
})

