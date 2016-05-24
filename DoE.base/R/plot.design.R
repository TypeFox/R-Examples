plot.design <- function(x, y=NULL, select=NULL, selprop=0.25, ask=NULL, ...){
      xnam <- deparse(substitute(x))
      if (!is.null(y)) ynam <- deparse(substitute(y))  ## for on the fly responses
      if (!"design" %in% class(x)){
         class(x) <- c("design", class(x)) 
         graphics::plot.design(x, y, ...)
         }
      else{## do the right thing for class design from package conf.design
         if (is.null(design.info(x)))
             graphics::plot.design(x, y,...)
         else{
             ## now designs generated with suite DoE.base etc.
             di <- design.info(x)
             nfac <- di$nfactors
             ## first process select
             ov <- FALSE ## overview by mosaic plots, ignoring response values
             if (is.null(select)) select <- names(di$factor.names)
             else if (is.list(select)){
                ov <- TRUE
                if (length(select)==1){
                    if (!all(select[[1]] %in% 1:nfac)) 
                       stop("length 1 list select must contain vector of factor position numbers")
                }
                else{
                ## all triples or quadruples from a subset of the factors
                if (select[[2]][1] %in% c("all2","all3","all4")){
                dimov <- as.numeric(substr(select[[2]],4,4))  ## pair, triple or quadruple
                select <- select[[1]]                         ## the selection of factors
                if (!is.numeric(select)) stop("if select is a list, its first element must be numeric")
                if (length(select)==1) {
                   hilf <- nchoosek(nfac-1, dimov-1)
                   hilf <- rbind(matrix(setdiff(1:nfac, select)[hilf], nrow=nrow(hilf)), select)
                   select <- as.list(as.data.frame(hilf))
                   }
                }
                else {
                    if (!all(sapply(select, function(obj) all(obj %in% 1:nfac))))
                        stop("select contains invalid entries")
                    if (!length(unique(sapply(select,length)))==1)
                        stop("all entries of list select must have the same length")
                    dimov <- length(select[[1]])
                }
                ov <- TRUE                                    ## overview mode (responses not used)
                }
             }
             else if (is.numeric(select)){ 
                if (!all(select %in% 1:di$nfactors))
                   stop("If numeric, select must contain integer numbers from 1 to the number of factors")
                select <- names(di$factor.names)[select]
             }
             if (is.character(select)){
                if (select[1] %in% c("all2","all3","all4","complete","worst","worst.rel","worst.parft","worst.parftdf")){
                  if (select[1] %in% c("all2","all3","all4")){
                    ## all triples or quadruples from all factors
                    ov <- TRUE
                    dimov <- as.numeric(substr(select,4,4))
                    select <- 1:di$nfactors
                  }
                  else{
                    ov <- TRUE
                    if (!select[1]=="complete"){
                      if (!is.numeric(selprop)){ 
                         warning("invalid selprop has been replaced by default 0.25")
                         selprop <- 0.25
                         }
                      if (!(selprop[1]>0 & selprop[1]<1)){
                         warning("invalid selprop has been replaced by default 0.25")
                         selprop <- 0.25
                      } 
                      selprop <- selprop[1]
                    }
                    type <- select[1]
                    select <- tupleSel(x, type=select[1], selprop=selprop) 
                        ## list of worst case tuples determined by tupleSel
                        if (length(select)==0){
                            if (type=="complete") stop("there are no tuples with complete confounding")
                            else stop("there are no tuples with confounding worse than the ", 1-selprop, " quantile")
                       }                 
                  }
                }
                else{
                if (all(select %in% Letters[1:di$nfactors]) & !all(names(di$factor.names) %in% Letters[1:di$nfactors]))
                    select <- names(di$factor.names)[which(Letters %in% select)]
                if (!all(select %in% names(di$factor.names)))
                    stop("select is invalid")
                }
             }
                         
             graphics <- FALSE
             table <- FALSE
             if (ov) {
                ## select contains the factor numbers to be used for the overview
                ## dimov contains the dimension of the plot 
                table <- TRUE
                if (is.list(select)){
                   tuples <- do.call(cbind, select)
                   select <- 1:nfac
                }
                else tuples <- nchoosek(length(select), dimov)
                askold <- devAskNewPage()
                if (is.null(ask)) ask <- dev.interactive(orNone=TRUE)
                devAskNewPage(ask=ask)
                suppressMessages(response.names(x) <- NULL)
                for (i in 1:ncol(tuples)) plot(x, select=select[tuples[,i]], ...)
                devAskNewPage(ask=askold)
                }
             else{
             if (is.null(di$quantitative)){
                 if(!(is.null(y) & is.null(di$response.names))) graphics <- TRUE
                 else table <- TRUE
             }
             else{ 
              if (all(is.na(di$quantitative)) | !any(di$quantitative==TRUE)){
                 if(!(is.null(y) & is.null(di$response.names))) graphics <- TRUE
                 else table <- TRUE
                 }
             }
             if (graphics){
                   if (is.null(y)) y <- x[,response.names(x)]
                   else {
                      if (is.character(y)){ 
                         if (!all(y %in% colnames(x))) stop("invalid names in y")
                         y <- x[,y]
                         }
                      }
                  if (is.data.frame(y)) y <- as.matrix(y)
                  if (!is.numeric(y)) stop("columns in y must be numeric")
                  ## as graphics::plot.design does not handle a matrix well
                  ## choose complicated way of handling responses 
                  if (!is.matrix(y)) {
                     y <- matrix(y, ncol=1)
                     if (!is.null(di$response.names)) 
                         colnames(y) <- di$response.names
                     else colnames(y) <- ynam
                   }
                   askold <- devAskNewPage()
                   if (is.null(ask) & ncol(y)>1) ask <- dev.interactive(orNone=TRUE)
                   devAskNewPage(ask=ask)
                   for (i in 1:ncol(y)){
                   cn <- colnames(y)
                   assign(cn[i], y[,i])
                   eval(parse(text=paste("graphics::plot.design(x[,c(select)],", 
                          cn[i], ", ask=ask, ...)")))
                   }
                   devAskNewPage(ask=askold)
                   }
             if (table){
                  ## process metric requests with special character strings 
                  ## from option sub
                  dots <- list(...)
                  if ("sub" %in% names(dots)){
                          if (dots[["sub"]] %in% c("GR","A", "rA", "sumPARFT", "sumPARFTdf")){
                          typ <- dots[["sub"]]
                          digits=4
                          if (typ=="GR") dots[["sub"]] <- paste("GR =",round(GR(x[,select])$GR, digits))
                          else {
                              l.s <- length(select)
                              if (l.s > 4) {
                                  dots[["sub"]] <- ""
                                  message("word length requests with option sub work for dimensions up to 4 only")
                                  }
                              else{ 
                              if (l.s==4){ 
                                 if (typ=="A"){
                                     l3 <- length3(x[,select])
                                     if (l3==0) dots[["sub"]] <- paste("A4 = ", round(length4(x[,select]), digits), sep="")
                                     else dots[["sub"]] <- paste("A3 = ", round(l3, digits), ", A4 = ", 
                                         round(length4(x[,select]), digits), sep="")
                                     }
                                 if (typ=="rA"){
                                     l3 <- length3(x[,select],rela=TRUE)
                                     if (l3==0) dots[["sub"]] <- paste("rA4 = ", round(length4(x[,select], rela=TRUE), digits), sep="")
                                     else
                                     dots[["sub"]] <- paste("rA3 = ", round(l3, digits), ", A4 = ", round(length4(x[,select]), digits), sep="")
                                 }
                                 if (typ=="sumPARFT"){
                                     l3 <- length3(x[,select])
                                     if (l3==0) dots[["sub"]] <- paste("sum(PARFT4) = ", round(attr(P4.4(x[,select],parft=TRUE),"sumPARFT4"),digits), sep="")
                                     else
                                     dots[["sub"]] <- paste("sum(PARFT3) = ", round(attr(P3.3(x[,select],parft=TRUE),"sumPARFT3"),digits), sep="")
                                 }
                                 if (typ=="sumPARFTdf"){
                                     l3 <- length3(x[,select])
                                     if (l3==0) dots[["sub"]] <- paste("sum(PARFTdf4) = ", round(attr(P4.4(x[,select],parftdf=TRUE),"sumPARFTdf4"),digits), sep="")
                                     else
                                     dots[["sub"]] <- paste("sum(PARFTdf3) = ", round(attr(P3.3(x[,select],parftdf=TRUE),"sumPARFTdf3"),digits), sep="")
                                 }
                                 }
                              if (l.s==3){ 
                                 if (typ=="A") dots[["sub"]] <- paste("A3 =", round(length3(x[,select]), digits))
                                 if (typ=="rA") dots[["sub"]] <- paste("rA3 =", round(length3(x[,select], rela=TRUE), digits))
                                 if (typ=="sumPARFT") 
                                    dots[["sub"]] <- paste("sum(PARFT3) = ", round(attr(P3.3(x[,select],parft=TRUE),"sumPARFT3"),digits), sep="")
                                 if (typ=="sumPARFTdf")
                                     dots[["sub"]] <- paste("sum(PARFTdf3) = ", round(attr(P3.3(x[,select],parftdf=TRUE),"sumPARFTdf3"),digits), sep="")
                                 }
                              if (l.s==2){ 
                                 if (typ=="A") dots[["sub"]] <- paste("A2 =", round(length2(x[,select]),digits))
                                 else{ dots[["sub"]] <- ""
                                 message("Relative words or projection average R^2s of length 2 have not been implemented.")}
                                 }
                          }
                      }
                   }
                   }
                   do.call(mosaic,c(list(table(x[,select])), dots))
             }
             if (!(table | graphics)) {if (is.null(y)) plot(undesign(x)[,select], ...)
                  else if (is.character(y)) plot(x[,c(select,y)], ...)
                       else plot(cbind(x[,select],y), ...)
             }
         }
      }
      }
}

