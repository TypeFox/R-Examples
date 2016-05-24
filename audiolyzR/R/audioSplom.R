audioSplom <- function (x = NULL, data, purge.plots = FALSE,
                        bins = 30, aspect = 1, radius = floor (sqrt(bins)) - 1,
                        key = "C", quality = "Major", tempo = 115,
                        directory = file.path (Sys.getenv("R_LIBS_USER"), "audiolyzR"),
                        output = file.path (tempdir(), "json_matrix"), write.to.home = NULL, ...)
{
  if (bins > 32) cat("WARNING: audiolyzR only supports 32 distinct notes, so /n bin values over 32 will generate strange results.")
  #checks whether audiolyzR is installed correctly
  audiolyzRcheck(directory)
  
  if( is.null(get("write.to.home", envir=asNamespace("audiolyzR"), inherits = FALSE) ) ) {
    message("R would like permission to write audioplot data files to: \n    ",
            file.path(Sys.getenv("HOME")),"  during this session.", 
            "\nOtherwise you will have to manually drag that folder to the synthesizer.",
            "\nIs that OK?")
    ANSWER <- readline("Type y or n")
    if (ANSWER %in% c("y","Y","yes","Yes","YES") )
      write.to.home <- TRUE else write.to.home <- FALSE
    unlockBinding("write.to.home", env = asNamespace("audiolyzR") )
    assign("write.to.home", write.to.home, envir = asNamespace("audiolyzR"))
  } else write.to.home <- get("write.to.home", envir = asNamespace("audiolyzR"), inherits=FALSE)
  
  #if a vector of variables is not provided
  if(is.null(x)){
    #obtain class of variables
    class.data <- data.frame(class=sapply(data,class), 
                             n.unique=sapply(data,function(x) length(unique(x)))) 
    #requires that at least 2 numeric variables exist
    if(nrow(class.data) < 3)
      stop("audioSplom requires at least 3 numeric variables")
    #take only numeric with more than 10 uniques
    data <- data[,which(class.data$class %in% c("numeric","integer") & class.data$n.unique >= 10)] 
    if(ncol(data) >= 20)
      stop("audioSplom only supports up to 20 variables")
  } 
  #otherwise use the given variables
  if(is.character(x)){
    if(length(x) < 3) stop("audioSplom requires at least 3 numeric variables")
    else data <- data[x] 
  }
  
  if(class(x)=="formula"){
    data <- as.data.frame(sapply(all.vars(x), function(y) data[y]))
  }
  
  print(hexplom(~data,xbins=bins,colramp = BTY,upper.panel = panel.hexboxplot, ...))
  
  mat.index <- combn(1:ncol(data),2) #combinations of plots for splom
  
  #defines reverb for each combination according to correlation
  reverb <- apply(mat.index[,(1:ncol(mat.index))],2,
                  function(x) abs(cor(data[x[1]],
                                      data[x[2]], use="complete.obs")))
  
  #applies the matricize function to each combination and joins the results into one data.frame
  matrices <- list()
  matrices <- apply(mat.index[,(1:ncol(mat.index))],2,
                    function(x) matrices[[which(mat.index[1,]==x[1] & mat.index[2,]==x[2])]] <- 
                      matricize(data[x[2]],data[x[1]],bins,radius,aspect,type="splom", ...))
  
  #remove rows that are all 0s
  matrices <- lapply(matrices, function(x) x[-which(x[3]==0),])  
  
  #extract label names for axes in the synth
  labels <- lapply(matrices, function(x) data.frame(main.x=unlist(strsplit(names(x[3]), "\\_"))[1],
                                                    main.y=unlist(strsplit(names(x[3]), "\\_"))[2]))   
  #store actual names in separate list
  for (i in 1:length(matrices))
    names(matrices[[i]]) <- c("Xcoord","Ycoord","plot.1","dist.1")
  
  global.options <- lapply(1:length(matrices), function(i)
  data.frame(key = paste(key), 
             quality = paste(quality), 
             reverb = reverb[i], 
             tempo = tempo, 
             maxdist = max(matrices[[i]][,4]),
             maxcount = max(matrices[[i]][,3])))
  
  #split each row frame into a list for each row
  matrices.main <- lapply(matrices, function(x) split(x, 1:nrow(x))) 
  
  #define overall dimensions of the json list
  if(aspect != 1) dimensions <- data.frame(X=max(matrices[[1]][1])+1,Y=max(matrices[[1]][2])+1) else
    dimensions <- data.frame(X=max(matrices[[1]][c(1,2)])+1,Y=max(matrices[[1]][c(1,2)])+1)
  
  #populate a list with all info for each graph
  max.out <- lapply(1:length(matrices), function(i)
    list(matrices.main[[i]],dimensions,global.options[[i]],labels[[i]]))
  
  for (i in 1:length(matrices)) {
    names(max.out[[i]]) <- c("matrix","dimensions","globaloptions","labels")
  }
    
  if (write.to.home)
    output <- file.path(Sys.getenv("HOME"),"json_matrix")
  
  if (purge.plots)
    unlink(output, recursive=TRUE)
  
  dir.create(output, showWarnings=FALSE)
  file <- lapply(labels, function(x) paste(x[1,2],"~",x[1,1],".json",sep=""))
  max.out <- lapply(max.out, function(x) toJSON(x, .withNames=TRUE))
  
  #write one file per plot to the json_matrix folder
  sapply(1:length(matrices), function(i) writeLines(max.out[[i]], paste(output,file[[i]],sep="/")))
  
  #launch synth
  systemRun(directory, output, write.to.home)
}
