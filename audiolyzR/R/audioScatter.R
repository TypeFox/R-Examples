audioScatter <-
function(x, y = NULL, z = NULL, data, purge.plots = FALSE, show.plots = TRUE,
         bins = 30, aspect = 1, radius = floor(sqrt(bins))-1, 
         key = "C", quality = "Major", tempo = 115,
         directory = file.path (Sys.getenv("R_LIBS_USER"), "audiolyzR"),
         output = file.path (tempdir(), "json_matrix"), write.to.home = NULL, ...)
  {
  if (bins > 32) 
    warning("audiolyzR only supports 32 distinct notes, so values over 32 may generate strange results.")
  
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
  
  #checks whether audiolyzR is installed correctly
  audiolyzRcheck(directory)
    
  #deparses the formula, if a formula is specified
  if(class(x)=="formula"){
    if(length(all.vars(formula)) > 3){
      stop("audioScatter only supports up to 3 variables. Try audioSplom instead.")
    } else {
      vars <- all.vars(x)
      x <- vars[2]
      y <- vars[1]
      if(length(vars) == 3)
        z <- vars[3]
      
    }
    
  }
  
  #Plots
  if(show.plots)
    if(!is.null(z)){
      print(hexbinplot(data[ ,y] ~ data[ ,x], xbins=bins, aspect=aspect, xlab=x, ylab=y, ...), position=c(  0, 0, .52, 1), more=TRUE)
      print(hexbinplot(data[ ,z] ~ data[ ,x], xbins=bins, aspect=aspect, xlab=x, ylab=z, ...), position=c( .48, 0, 1, 1))
    } else {
      print(hexbinplot(data[ ,y] ~ data[ ,x], xbins=bins, aspect=aspect, xlab=x, ylab=y, ...))
    }
  
  #global synth variables for Audiolize
  if(!is.null(z)) {
    reverb <- abs(cor(data[,z],data[,x], use="complete.obs"))} else {
      reverb <- abs(cor(data[,y],data[,x], use="complete.obs"))}
  
  #creates matricized frames for ea var
  mat <- matricize(data[y],data[x], bins, radius, aspect, type = "scatter", ...) 
  
  #global dimensions for Audiolize
  if(aspect != 1) dimensions <- data.frame(X=max(mat[,1])+1,Y=max(mat[,2])+1) else
    dimensions <- data.frame(X=max(mat[,c(1,2)])+1,Y=max(mat[,c(1,2)])+1)
  global.options <- data.frame(key=paste(key), 
                               quality=paste(quality), 
                               reverb=reverb, 
                               tempo=tempo,
                               maxdist=max(mat[,4]),
                               maxcount=max(mat[,3]))
  
  #performs the same matricize process on x if a third variable is specified, joins with mat
  if(!is.null(z)){
    mat.dual <- data.frame(mat, matricize(data[z],data[x], bins, radius, aspect, type = "scatter", ...)[3:4])
    mat.dual <- mat.dual[-which(mat.dual[,3]==0 & mat.dual[,5]==0),] #gets rid of zero components
    #stores actual names
    labels <- data.frame(main.x=paste(x),main.y=paste(y),cond.y=paste(z)) 
    #rgeneric names
    names(mat.dual) <- c("Xcoord","Ycoord","plot.1","dist.1","plot.2","dist.2") 
    #turns into a list by row for RJSON export
    mat.dual.list <- split(mat.dual, 1:nrow(mat.dual))
    max.out <- list(mat.dual.list,dimensions,global.options,labels)
    names(max.out) <- c("matrix","dimensions","globaloptions","labels")
  } else {       
    mat <- mat[-which(mat[,3]==0),]
    labels <- data.frame(main.x=x,main.y=y) #stores actual names
    names(mat) <- c("Xcoord","Ycoord","plot.1","dist.1") #replaces with generic names
    mat.main <- split(mat, 1:nrow(mat))
    max.out <- list(mat.main,dimensions,global.options,labels)
    names(max.out) <- c("matrix","dimensions","globaloptions","labels")
  }
  
  if (write.to.home)
    output <- file.path(Sys.getenv("HOME"),"json_matrix")
  
  if (purge.plots)
    unlink(output, recursive=TRUE)
  
  dir.create(output, showWarnings=FALSE)
  if(!is.null(z)) file <- paste(labels[1,1],"~",labels[1,2],"+",labels[1,3],".json",sep="")
  if(is.null(z)) file <- paste(labels[1,1],"~",labels[1,2],".json",sep="")
  max.out <- toJSON(max.out, .withNames=TRUE)
  writeLines(max.out, file.path(output,file) )
  
  systemRun(directory, output, write.to.home)
 
}
