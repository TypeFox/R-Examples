# # logLik function:
# logLik.mlVAR <- function(object){
#   res <- object$pseudologlik
#   class(res) <- "logLik"
#   attr(res, "df") <- object$df
#   return(res)
# } 

# Plot function:
plot.mlVAR <- function(x, type = c("fixed","se","random","subject"), lag = 1,subject,...){
  if (type[[1]]=="subject" & missing(subject)){
    stop("'subject' is needed to plot individual network")
  }
  
  Nodes <- rownames(x$fixedEffects)
  if (type[[1]]=="fixed"){
    
    Net <- x$fixedEffects[,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(x$fixedEffects))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
  } else if (type[[1]]=="se"){
    
    Net <- x$se.fixedEffects[,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(x$se.fixedEffects))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
  } else if (type[[1]]=="random"){
    
    Net <- x$randomEffectsVariance[,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(x$randomEffectsVariance))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
    
  }  else if (type[[1]]=="subject"){
    
    Net <- x$randomEffects[[subject]][,grepl(paste0("^L",lag,"_(",paste(Nodes,collapse="|"),")$"),colnames(x$randomEffects[[subject]]))]
    Graph <- qgraph(t(Net), labels=rownames(Net), ...)
    
  } else stop("'type' is not supported")
  
  invisible(Graph)
}


# Print and summary:
summary.mlVAR <- function(object,...){
  input <- object$input
  
  inputstr <- paste(sapply(seq_along(input),function(i)paste0(names(input)[i],":\t\t",paste(input[[i]],collapse=", "))), collapse = "\n")
  
  cat(paste0("==== mlVAR results ====\n",inputstr,"\n\nNumber of random effects:\t\t",length(object$randomEffects),"\n",
             "pseudo log-likeligood:\t\t",round(object$pseudologlik,2),"\n",
             "Degrees of Freedom:\t\t",round(object$df,2),"\n",
             "BIC:\t\t",round(object$BIC,2)             
             ))
}

print.mlVAR <- function(x,...) summary.mlVAR(x,...)
