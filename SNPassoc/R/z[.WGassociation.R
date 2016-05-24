`[.WGassociation` <-
function(x,i,j,  ...){
    if (missing(i)) i<-1:nrow(x)
    if (missing(j)) j<-1:ncol(x)

    if (is.numeric(i)) i<-(1:nrow(x))[i]
    if (is.numeric(j)) j<-(1:ncol(x))[j]

    if (is.character(i)) i<-match(i,rownames(x))
    if (is.character(j)) j<-match(j,colnames(x))
    
    if (is.logical(i)) i<-(1:nrow(x))[i]
    if (is.logical(j)) j<-(1:ncol(x))[j]


    if (!(1 %in% j)) j<-c(1,j)
    if (length(j)==1) if(j==1) j<-1:2
    if (any(!(i %in% 1:nrow(x)))) stop("Undefined rows selected")
    if (any(!(j %in% 1:ncol(x)))) stop("Undefined cols selected")
    
    out<-attr(x,"pvalues")[i,j]  # data.frame

    attr(out, "tables")         <- attr(x, "tables")[i]
    attr(out, "label.SNPs")     <- attr(x, "label.SNPs")[i]
    attr(out, "colSNPs")        <- attr(x, "colSNPs")[i]
    attr(out, "gen.info")       <- attr(x, "gen.info")[i,]
    attr(out, "whole")          <- attr(x, "whole")
    attr(out, "pvalPerm")       <- attr(x, "pvalPerm")[i,]
    attr(out, "pvalues")        <- out
    attr(out, "models")         <- c(0,attr(x, "models"))[j]
    attr(out, "quantitative")   <- attr(x, "quantitative")
    attr(out, "fast")           <- attr(x, "fast")
    class(out) <-  c("WGassociation", "data.frame")

    out }

