#' Transforms variables to cases
#'
#' Creates a transfromed xpssFrame.
#'
#' @usage xpssVarsToCases(x, from, idVar = NULL, indexVar = NULL, nullArg = "keep", 
#' countVar = NULL, varLabels = list(id = NULL, index = NULL, count = NULL))
#' @param x as a (non-empty) data.frame or input data of class "xpssFrame". 
#' @param from variable that opens the span.
#' @param idVar determines whether an id-variable should be created.
#' @param indexVar determines whether an index-variable should be created.
#' @param nullArg Can be either "keep" or "drop".
#' @param countVar determines whether a counter should be created?
#' @param varLabels determines whether labels for id-, index- and count variables are set.
#' @return Returns the transformed xpssFrame.
#' @author Andreas Wygrabek
#' @importFrom plyr ldply llply
#' @importFrom tidyr gather
#' @examples 
#' data(fromXPSS)
#'
#' xpssVarsToCases(fromXPSS, from = list(c("newVar", "V7_1, V7_2")), 
#' idVar = "myID", indexVar = "myIndex", nullArg = "drop", countVar = "Counter")
#' @export
xpssVarsToCases <- function(x, from, idVar = NULL, indexVar = NULL, nullArg = "keep", countVar = NULL, varLabels = list(id = NULL, 
                                                                                                                      index = NULL,
                                                                                                                        count = NULL)){
options(warn = -1)

####################################################################
####################################################################

functiontype <- "DM"
x <- applyMetaCheck(x)

####################################################################
####################################################################
####################################################################


if(length(from) < 1){
    stop("Argument from is required")
}

testLengthFrom <- sapply(from, function(x){
    length(x)
})

if(sum(testLengthFrom)/length(from) != 2){
    stop("Each element of from has to be length 2")
}  


if(!is.null(varLabels[["id"]]) & is.null(idVar)){
    stop("Create an idVar first, before you label it")
}

if(!is.null(varLabels[["index"]]) & is.null(indexVar)){
    stop("Create an indexVar first, before you label it")
}

if(!is.null(varLabels[["count"]]) & is.null(countVar)){
    stop("Create a countVar first, before you label it")
}

vars <- llply(.data = from, 
              function(x){
                  x[[2]]
              })

x[,"idVar"] <- 1:nrow(x)

varnames <- llply(.data = from, 
                  function(x){
                      x[[1]]
                  })

eval(parse(text = paste("ini <- gather(x,'key','make'",",",unlist(vars)[1],")", sep = "")))
ini <- ini[,-which(colnames(ini)=="key" | colnames(ini)=="make")]

eval(parse(text = paste("temp",1:length(from)," <- gather(x,'key','make'",",",unlist(vars),")[,'make']", sep = "")))
mat <- do.call("cbind", as.list(parse(text = paste0("temp",1:length(vars)))))

colnames(mat) <- varnames

out <- cbind(ini, mat)

if(nullArg == "drop"){
    if(length(colnames(mat)) > 1){
    ind <- rowSums(mat, na.rm = TRUE)
    ind <- which(!is.na(ind))} else {
        ind <- which(!is.na(mat[,1]))
    }
} else if(nullArg == "keep") {
    ind <- rep(TRUE,nrow(out))
}

if(!is.null(countVar) & nullArg == "keep"){
    
    counts <- ldply(split(out, out[,"idVar"]),
                    function(x){
                        nrow(x)
                    })
    
     out[,countVar] <- apply(out,1,function(x){
         
         POS <- which(as.integer(x["idVar"]) == 1:nrow(counts)) 
         counts[POS,2]
     })
} else if(!is.null(countVar) & nullArg == "drop") {
    
    counts <- ldply(split(out[ind,], out[ind,"idVar"]),
                    function(x){
                        nrow(x)
                    })
    
    out[,countVar] <- apply(out,1,function(x){
        
        POS <- which(as.integer(x["idVar"]) == 1:nrow(counts)) 
        counts[POS,2]
    })
}

if(!is.null(idVar)){
    out <- out[,c("idVar", colnames(out)[-length(colnames(ini))])]
    out <- xpssSortCases(out,"idVar","A")
    colnames(out)[colnames(out) == "idVar"] <- idVar 
} else {
    out <- xpssSortCases(out,"idVar","A")
    out[,which(colnames(out) == "idVar")] <- NULL
}

if(!is.null(indexVar)){
    out[,indexVar] <- rep(1:(nrow(out)/nrow(x)),length(unique(out[,1])))
}

attrDF <- attributes(out)
attrVars <- sapply(out, function(x){
    attributes(x)
})
orderRN <- rownames(out[(rownames(out) %in% ind) | ind == TRUE,])

out <- out[(rownames(out) %in% ind) | ind == TRUE,]

attrDF$row.names <- orderRN
attributes(out) <- attrDF
eval(parse(text = paste("attributes(out[,",1:length(colnames(out)),"]) <- attrVars[[",1:length(colnames(out)),"]]", sep = "")))
rownames(out) <- orderRN

attr_backup <- attributesBackup(x)

pos <- which(!names(attr_backup$global) %in% names(attributes(out)))

for(i in 1:length(names(attr_backup$global)[pos])) {
  attr(out,names(attr_backup$global)[pos][i]) <- attr_backup$global[pos][[i]]
}

class(out) <- c("xpssFrame", "data.frame")

if(!is.null(idVar) & !is.null(varLabels[["id"]])){
    out <- xpssVariableLabels(out, variables = c(idVar), labels = varLabels[[1]])}

if(!is.null(indexVar) & !is.null(varLabels[["index"]])){
    out <- xpssVariableLabels(out, variables = c(indexVar), labels = varLabels[[2]])}

if(!is.null(countVar) & !is.null(varLabels[["count"]])){
    out <- xpssVariableLabels(out, variables = c(countVar), labels = varLabels[[3]])}

options(warn = 0)

return(out)
}
