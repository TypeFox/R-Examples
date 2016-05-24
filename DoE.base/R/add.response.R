add.response <- function(design, response, rdapath=NULL, replace=FALSE, 
     InDec=options("OutDec")[[1]], tol=.Machine$double.eps ^ 0.5, ...){
     ## invalid paths are reasonably warned about outside of this program
     if (!is.null(rdapath)){
        load(rdapath)
        if (!exists(design)) stop(design, " not found in ", rdapath)
        ## desnam <- design  ## that is not needed, is it ?
        assign(design, get(design))  ## why this?
        design <- get(design)
        if (!"design" %in% class(design)) stop("design must be of class design")
        }
     if (is.character(response)){ 
        if (!substr(response, nchar(response)-2,nchar(response))=="csv") 
            stop("response must be a numeric vector or matrix, a data frame or the name of a csv-file")
        if (InDec==".") assign("response", read.csv(response))
             else assign("response", read.csv2(response))
     }

### take care of name that is automatically stored and should subsequently only be added if 
### it has been changed in the mean time

    rn <- make.names(deparse(substitute(response)))    ## make.names makes valid R names from e.g. rnorm(8)

    if (!"design" %in% class(design)) stop("add.response works on class design objects only.")
    di <- design.info(design)
    
    if (!(is.numeric(response) | is.data.frame(response)))
        stop("response must be a numeric vector or matrix or a data frame.")

    ## response to become a data frame with reasonable column names
    if (is.matrix(response)) {
        if (ncol(response) == 1 & is.null(colnames(response))) colnames(response) <- rn
        response <- data.frame(response)
        }
    if (is.numeric(response)) {
            response <- data.frame(response)
            colnames(response) <- rn
            }

    if (is.data.frame(response)){
        numc <- sapply(response, is.numeric)
        if (!any(numc))
            stop("response does not contain any numeric variables.")
            }
        else if (is.matrix(response)) numc <- rep(TRUE,ncol(response))
        else numc <- TRUE
    if (is.data.frame(response) | is.matrix(response))
        if ("run.no" %in% colnames(response)){ 
            rrno <- response[,"run.no"]
            drno <- run.order(design)$run.no
            if (is.factor(drno)) response[,"run.no"] <- rrno <- factor(rrno, levels=levels(drno))
            else if (is.factor(rrno)) response[,"run.no"] <- rrno <- as.numeric(as.character(rrno))
            if (!all(rrno==drno)){
               design <- design[ord(matrix(as.numeric(run.order(design)$run.no),ncol=1)),]
               response <- response[ord(matrix(as.numeric(response[,"run.no"]),ncol=1)),]
            }}
    
    fnam <- names(di$factor.names)
    blocknam <- di$block.name
    if (!nrow(response) ==nrow(design))
            stop("wrong number of observations in response")
    respnam <- setdiff(colnames(response)[which(numc)], c(fnam, blocknam, colnames(run.order(design)), "Name", "name"))
    if (length(respnam)==0) stop("no response variables in response")
    
     ## check consistency in case of overlapping variables
     ## NA variables must be removed from this, in order to allow adding only some response variables comfortably
     ## response variables are checked later
       ## this part also prevents replacement for ill-sorted designs, if the Factor variables are also kept
       
       ## gleich contains the non-NA variables occurring in both design and hilf but not in respnam
     gleich <- setdiff(intersect(colnames(design)[which(!sapply(design, function(.x) all(is.na(.x))))],colnames(response)),respnam)
     if (length(gleich)>0){
       if (!all(response[,gleich]==design[,gleich])){
            fn <- intersect(fnam,gleich) 
            if (length(fn) > 0) 
                for (ffn in fn) if (all(as.character(design[[ffn]])==as.character(response[[ffn]]))){ 
                          response[[ffn]] <- design[[ffn]]
                          gleich <- setdiff(gleich, ffn)
                          }
                          else{ 
                              if (isTRUE(all.equal(as.numeric(as.character(design[[ffn]])),as.numeric(as.character(response[[ffn]])),tolerance=tol))){ 
                                response[[ffn]] <- design[[ffn]]
                                gleich <- setdiff(gleich, ffn)
                                }
                             else{
                                print(ffn)
                                print(cbind(as.character(design[[ffn]]), as.character(response[[ffn]])))
                                stop("There are variables of same names but with different content in design and response")
                             }
                          }
            }
            }
    
    ## treat response variables
    ## if overlap and design contains any response neither equal to update nor NA --> prevent replacement
       ## now gleich contains the new names that are already in the design
    gleich <- intersect(colnames(design),respnam)
    if (length(gleich) > 0 & !replace){
        ## check non-NA values for approximate equality
              dm <- as.matrix(design[,gleich])
              rm <- as.matrix(response[,gleich])
        if (!(all(is.na(dm)) | 
               isTRUE(all.equal(dm[!is.na(dm)],rm[!is.na(dm)],tolerance=tol))))
           stop("response variables ", paste(gleich, sep=", "), " already exist in the design with partly different values!")
        if (!identical(dm[!is.na(dm)],rm[!is.na(dm)])){
              ## replace rm values with small numeric differences (all.equal was TRUE) by dm values
              ## necessary because of inadequate forced rounding in software like Excel 
              ##             for calculated responses with many decimal places
              rm[!is.na(dm)] <- dm[!is.na(dm)]
              response[,gleich] <- as.data.frame(rm)
            }
        }
    design[,respnam] <- response[,respnam]
    attr(design,"desnum") <- cbind(desnum(design), as.matrix(response[,respnam,drop=FALSE]))
    attr(design,"design.info")$response.names <- union(response.names(design), respnam)
         ## union for the case that some added response overwrite existing ones
         ## if no existing ones, NULL causes no problem
    design
}