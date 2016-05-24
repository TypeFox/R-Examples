# conversion to class 'cross'
# return F2 intercross

gpData2cross <- function(gpData,...){
     # check for class
     if(class(gpData)!="gpData") stop("object '",substitute(gpData),"' not of class 'gpData'")  
     # check on geno and map
     if(is.null(gpData$geno) | is.null(gpData$map)) stop("'geno' and 'map' needed in",substitute(gpData))
     if(dim(gpData$pheno)[3] > 1) stop("You can only use unreplicated values for cross!")
     else{
       # use codeGeno if not yet done
       if(!gpData$info$codeGeno) stop("Use function codeGeno before gpData2cross!")
     
       # only use individuals with genotypes and phenotypes
       genoPheno <- gpData$covar$id[gpData$covar$genotyped & gpData$covar$phenotyped]
       # read information from gpData
       geno <- data.frame(gpData$geno[rownames(gpData$geno) %in% genoPheno,])
       phenoDim <- dim(gpData$pheno)
       phenoNames <- dimnames(gpData$pheno)
       phenoDim[1] <- sum(dimnames(gpData$pheno)[[1]] %in% genoPheno)
       phenoNames[[1]] <- dimnames(gpData$pheno)[[1]][dimnames(gpData$pheno)[[1]] %in% genoPheno]
       pheno <- gpData$pheno[dimnames(gpData$pheno)[[1]] %in% genoPheno, ,]
       pheno <- array(pheno, dim=phenoDim)
       dimnames(pheno) <- phenoNames
       pheno <- apply(pheno, 2, rbind) # possible because of unreplicated data!!!
       rownames(pheno) <- phenoNames[[1]]
       if(dim(gpData$pheno)[3]>1) pheno$repl <- rep(1:dim(gpData$pheno)[3], each=dim(gpData$pheno)[1])
       if(!is.null(gpData$phenoCovars)) pheno <- cbind(pheno, data.frame(apply(gpData$phenoCovars[dimnames(gpData$phenoCovars)[[1]] %in% genoPheno, ,], 2, rbind)))
       map  <- gpData$map
       n <- nrow(geno)
       pheno <- as.data.frame(pheno)
     }
     
     # split markers (+pos) and genotypes on chromosomes
     genoList <- split(cbind(rownames(map),map$pos,t(geno)),map$chr)
     # result is a list       
     # function to bring each list element in right format
     addData <- function(x){                       
         ret <- list()
         Nm <- length(x)/(n+2)
         
         # elements of x:
         #  1:Nm: marker names
         #  (Nm+1):(2*Nm): marker positions
         #  rest: genotypes as vector
         
         # add 1 to genotypes
         # coding for F2 intercross: AA=1, AB=2, BB=3
         ret[["data"]] <- matrix(as.numeric(x[-(1:(2*Nm))])+1,nrow=n,ncol=Nm,byrow=TRUE,dimnames=list(NULL,x[1:Nm]))
         ret[["map"]]  <- as.numeric(x[(Nm+1):(2*Nm)])
         names(ret[["map"]]) <- x[1:Nm]
         # this may have to be changed
         class(ret) <- "A"
         ret
     }
     #                  
     # apply function to each list element
     genoList <- lapply(genoList,addData)
     # create object 'cross'
     cross <- list(geno=genoList,pheno=pheno)
     class(cross) <- c("f2","cross")
     cross                   
}
