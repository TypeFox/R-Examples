print.taglist <- function(x, LineLength = getOption("width"), offset.start = 0,
                          withFinalLineBreak = TRUE, first.print = NULL,
                          ErrorOrWarn = "warn", ...){
   xc <- as.character(deparse(substitute(x)))
   ll <- length(x)
   LineL <- max(LineLength-2,0)
   LineBreak <- NULL
   mi50 <- min(LineLength,50)
   maL <- max(3*LineLength,getOption("width"))
   if(ll){
      offS <- paste(rep(" ", offset.start), collapse = "")
      for(i in 1:ll){
          trystr0 <- paste(names(x[i]),x[[i]],sep = "=")
          if(i==1){
            actstr <-  trystr  <- trystr0
            trystr0 <- NULL
            if(length(first.print))
               cat(first.print)
          }else{
             trystr  <- paste(actstr, trystr0, sep = ",")
          }

          ntry <- nchar(trystr) + offset.start
          if (ntry < LineL){
              actstr <- trystr
          }else{
              WarnOrErrorFct <- if(pmatch(tolower(ErrorOrWarn),
                                   c("warn","error"))==1)
                                   warning else stop
              if(ntry > maL) WarnOrErrorFct(gettextf(
                      "Some elements of %s are too long",
                       if(nchar(xc)> mi50) paste(substr(xc,1,mi50),"...") else
                                xc))
              if(actstr!=offS) cat(LineBreak, actstr,",%",sep="")
              LineBreak <- "\n"
              actstr <- paste(offS, trystr0, sep = "")
          }
      }
      if(nzchar(actstr)) {
         if(i>1) cat("\n")
         cat(actstr, sep="")
         if(withFinalLineBreak) cat("\n")
      }
   }
   return(invisible())
}

taglist <- function(..., list = NULL, defname = "V"){
  dots <-  c(list,match.call(call = sys.call(),
                       expand.dots = FALSE)$"...")
  ldots <- length(dots)
  if(ldots){
     defname <- unique(defname)
     defnames <- if(length(defname)<ldots){
            paste(rep(defname,length.out=ldots), 1:ldots, sep = "")} else defname
     nms <- names(dots)
     if(is.null(nms)) nms <- rep("", ldots)
     nms[ nms == "" ] <- defnames[nms == ""]
     names(dots) <- nms
     return(structure(as.list(dots), class = c("taglist","list")))
  }
  return(invisible())
}
