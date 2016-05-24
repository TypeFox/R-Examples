`snp` <-
function (x, sep = "/", name.genotypes, reorder="common", remove.spaces = TRUE, 
       allow.partial.missing = FALSE) 
{

if (is.snp(x))
 {
  object<-x
 }

else
 { 
  if (sum(is.na(x)) == length(x))
     {
       object<-rep(NA, length(x))      
       attr(object, "allele.names") <- NULL
       class(object) <- c("snp","logical")
       return(object)

     }

 if(missing(name.genotypes))
 {
    alleles<-NULL
    x.d <- dim(x)
    x <- as.character(x)
    dim(x) <- x.d
    x[is.na(x)] <- ""
   
    if (remove.spaces) {
        xdim <- dim(x)
        x <- gsub("[ \t]", "", x)
        dim(x) <- xdim
        
    }


    if (!is.null(dim(x)) && ncol(x) > 1) 
        parts <- x[, 1:2]
    else {
        if (sep == "") 
            sep <- 1
        if (is.character(sep)) {
            part.list <- strsplit(x, sep)
            part.list[sapply(part.list, length) == 0] <- NA
            half.empties <- lapply(part.list, length) == 1
            part.list[half.empties] <- lapply(part.list[half.empties], 
                c, NA)
            empties <- is.na(x) | lapply(part.list, length) == 
                0
            part.list[empties] <- list(c(NA, NA))
            parts <- matrix(unlist(part.list), ncol = 2, byrow = TRUE)
        }
        else if (is.numeric(sep)) 
            parts <- cbind(substring(x, 1, sep), substring(x, 
                sep + 1, 9999))
        else stop(paste("I don't know how to handle sep=", sep))
    }
    mode(parts) <- "character"
    temp <- grep("^[ \t]*$", parts)
    parts[temp] <- NA
    if (!allow.partial.missing) 
        parts[is.na(parts[, 1]) | is.na(parts[, 2]), ] <- c(NA, 
            NA)
    alleles <- unique(c(na.omit(parts)))

    if(length(alleles)>2)
       stop("SNP must have only two alleles")


    tmp <- ifelse(is.na(parts[, 1]) & is.na(parts[, 2]), NA, 
        apply(parts, 1, paste, collapse = "/"))
    object <- factor(tmp)
    
    ll <- levels(object) <- na.omit(levels(object))

    if (length(ll)==4)
     {
       object[object==ll[3]]<-ll[2]
       object<-factor(object) 
     }
   
    control <- paste(rep(alleles[1],2),collapse="/")%in%ll
 
    if (sum(control)==0 & length(ll)==3)
     {
       object[object==ll[2]]<-ll[1]
       object<-factor(object) 
     }


    control <- paste(rep(alleles[2],2),collapse="/")%in%ll
 
    if (sum(control)==0 & length(ll)==3)
     {
       object[object==ll[3]]<-ll[2]
       object<-factor(object) 
     }

    if (length(object)==sum(is.na(object)))
     stop("choose the correct character separator to divide alleles") 

    class(object) <- c("snp","factor")
    object<-reorder.snp(object, ref=reorder)
    attr(object, "allele.names") <- alleles
}

else
 {
   if (any(is.na(match(x[!is.na(x)],name.genotypes))))
    stop("'name.genotypes' must match with the observed genotypes")
   x[x==name.genotypes[1]]<-"A/A"
   x[x==name.genotypes[2]]<-"A/B"
   x[x==name.genotypes[3]]<-"B/B"
   object<-as.factor(x)
   attr(object, "allele.names") <- c("A","B")
   class(object) <- c("snp","factor")
  }
 }

object

}

