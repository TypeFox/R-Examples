#' Write formatted data file
#' 
#' Write data contained in a list of dataframes in a single file (NONMEM/Monolix format)
#' or in several files as tables
#' 
#' See http://simulx.webpopix.org/mlxr/writedatamlx/ for more details.
#' @param r a list of dataframes
#' @param result.file a string with the name of the file 
#' @param result.folder a string with the name of the folder 
#' @param sep (default = ",") 
#' @param ext a string with the extension of the file names 
#' @param digits (default = 5) 
#' @param app.file  TRUE/FALSE (default=FALSE) append to file 
#' @param app.dir  TRUE/FALSE (default=FALSE) append to dir 
#' @examples
#' \dontrun{
#' modelPK <- inlineModel("
#' [LONGITUDINAL]
#' input = {V, Cl, a1}
#' EQUATION:
#' Cc = pkmodel(V, Cl)
#' DEFINITION:
#' y1 ={distribution=lognormal, prediction=Cc, sd=a1}
#' ")
#' adm  <- list(amount=100, time=seq(0,50,by=12))
#' p <- c(V=10, Cl=1, a1=0.1)
#' y1 <- list(name=c('y1'), time=seq(5,to=50,by=5))
#' res <- simulx(model=modelPK, treatment=adm, parameter=p, output=y1)
#' writeDatamlx(res, result.file="res.csv")
#' writeDatamlx(res, result.file="res.txt", sep="\t")
#' writeDatamlx(res, result.folder="res")
#' }
#' @importFrom utils write.table
#' @export
writeDatamlx <- function(r,result.file=NULL,result.folder=NULL,sep=",",ext=NULL,digits=5,app.file=F,app.dir=F) 
{
  if (!is.null(result.folder)){
    if (app.dir==F){
      unlink(result.folder, recursive = TRUE, force = TRUE)
      Sys.sleep(0.1)
      dir.create(result.folder, showWarnings = FALSE, recursive = FALSE, mode = "0777")
    }
    
    nr <- names(r$parameter)
    r0 <- r
    r0[nr] <- NULL
    r0["group"] <- NULL
    nr0 <- names(r0)
    f.names <- NULL
    if (is.null(ext)){
      if (sep %in% c(",",";"))
        ext <- ".csv"
      else
        ext <- ".txt"
    }
    for (k in (1:length(nr0))){
      rk <- r[[nr0[k]]]
      if (!is.null(ncol(rk))){
        fk <- file.path(result.folder,paste0(nr0[k],ext))
        for (j in (1:ncol(rk))){
          if (typeof(rk[,j])=="double")
            rk[,j] <- round(rk[,j], digits=digits)
          i.na <- which(is.na(rk[,j]))
          if (length(i.na)>0)
            rk[i.na,j]="."
        }
        
        if (app.file == F) 
          write.table(rk,fk,row.names=FALSE,quote=FALSE,sep=sep,append=F)
        else
          write.table(rk,fk,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=sep,append=T)
      }
    }
  }
  
  if (!is.null(result.file)){
    y.attr <- sapply(r,attr,"type")
    j.long <- which(y.attr=="longitudinal")
    y <- NULL
    #     if (length(j.long)==1){
    #       y <- r[[j.long]]
    #     }else if (length(j.long)>1){
    for (k in (1:length(j.long))){
      rk <- r[[j.long[k]]]
      nk <- names(r[j.long[k]])
      names(rk)[names(rk)==nk] <- "y"
      if (length(j.long)==1){
        y <- rk
      } else {
        yk <- cbind(rk,list(ytype=k))
        y <- rbind(y,yk)
      }
    }
    M <- y
    
    if (!is.null(r$treatment)){
      trt <- r$treatment
      trt$y <- NA
      M <- merge(M,trt,all=TRUE)
    }
    
    j.reg <- which(y.attr=="regressor")
    if (length(j.reg)>0){
      for (k in (1:length(j.reg)))
        M <- merge(M,r[[j.reg[k]]],all=TRUE)
    }
    
    
    if (!is.null(r$parameter)){
      n1 <- ncol(M)
      M <- merge(M,r$parameter,all=T)
      n2 <- ncol(M)
      occ <- M[,(n1+1):n2]
      n <- nrow(occ)
      for (i in (2:n)){
        if (any(is.na(occ[i,])))
          occ[i,] <- occ[(i-1),]
      }
      M[,(n1+1):n2] <- occ
    }
    
    for (k in (1:ncol(M))){
      if (typeof(M[,k])=="double")
        M[,k] <- round(M[,k], digits=digits)
      i.na <- which(is.na(M[,k]))
      if (length(i.na)>0)
        M[i.na,k]="."
    }
    
    lo <- NULL
    if (!is.null(M$pop))   lo <- c(lo, "pop")
    if (!is.null(M$rep)) lo <- c(lo, "rep")
    if (!is.null(M$id)) lo <- c(lo, "id")
    if (!is.null(M$time)) lo <- c(lo, "time")
    if (!is.null(M$ytype)) lo <- c(lo, "ytype")
    lo <- paste(lo,collapse=",")
    eval(parse(text=paste0("M <- M[with(M, order(",lo,")), ]")))
    
    if (app.file == F) 
      write.table(M,result.file,row.names=FALSE,quote=FALSE,sep=sep,append=F)
    else
      write.table(M,result.file,row.names=FALSE,col.names=FALSE,quote=FALSE,sep=sep,append=T)
    
    #     write.table(M,result.file,row.name=FALSE,quote=FALSE,sep=sep)
    #     zz <- file(result.file)
    #     close(zz)
  } 
}