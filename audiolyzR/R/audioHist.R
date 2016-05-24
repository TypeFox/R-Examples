audioHist <- function(x, name = "Variable", purge.plots = FALSE,
                      bins = 30, breaks = "Scott", radius = floor(sqrt(bins))-1, 
                      key = "C", quality = "Major", tempo = 80, reverb = 1,
                      directory = file.path (Sys.getenv("R_LIBS_USER"), "audiolyzR"),
                      output = file.path (tempdir(), "json_matrix"), write.to.home = NULL, ...) 
{
  if (bins > 32) cat("WARNING: audiolyzR supports 32 distinct notes, so /n bin values over 32 can generate hard-to-hear results.")
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
  
  h <- hist(x, breaks=breaks ,freq=TRUE, xlab=name, main=paste("Histogram of", name))
  bins <- ceiling(bins/length(h$counts)) * length(h$counts)
  counts <- ceiling(rescale(h$counts, c(0,bins-1)))
  
  mat <- matrix(sapply(counts, function(x)
    rep( c(rep(0, ( (bins) - x) ), rep(1, x) ), bins/length(counts))), ncol = bins)
  
  coords <- data.frame(Xcoord = rep(seq(1:bins)-1, bins),
                       Ycoord = rev(rep(c(0:(bins-1)), each=bins)))
  
  coords$V3 <- sapply(1:nrow(coords), function(x) 
    mat[coords[x ,2]+1, coords[x ,1]+1] )
  
  mat <- matricize (a = coords, b = name, bins = bins, 
                    radius = radius, type = "hist")
  
  dimensions <- data.frame(X = max(mat[,1])+1,
                           Y = max(mat[,2])+1)
  global.options <- data.frame(key=paste(key), 
                               quality=paste(quality), 
                               reverb=reverb, 
                               tempo=tempo,
                               maxdist=max(mat[,4]),
                               maxcount=max(mat[,3]))
  mat <- mat[-which(mat[,3]==0),]
  labels <- data.frame(main.x=paste(name),main.y="Freq") #stores actual names
  names(mat) <- c("Xcoord","Ycoord","plot.1","dist.1") #replaces with generic names
  mat.main <- split(mat, 1:nrow(mat))
  max.out <- list(mat.main,dimensions,global.options,labels)
  names(max.out) <- c("matrix","dimensions","globaloptions","labels")
  
  if (write.to.home)
    output <- file.path(Sys.getenv("HOME"),"json_matrix")
  
  if (purge.plots)
    unlink(output, recursive=TRUE)
  
  dir.create(output, showWarnings=FALSE)
  file <- paste(labels[1,1],"_",labels[1,2],".json",sep="")
  max.out <- toJSON(max.out, .withNames=TRUE)
  writeLines(max.out, paste(output,file,sep="/"))
  
  return(invisible(list(hist(x, breaks="Scott", freq=TRUE, 
                             xlab=name, main=paste("Histogram of", name)),
                        systemRun(directory, output, write.to.home))))
}
