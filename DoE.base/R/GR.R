GR <- function(ID, digits=2){
    ## function for calculating generalized resolution
    ## up to resolution V (5 means: 5plus)
    
    ## RPFTs are returned for up to resolution IV, otherwise NULL
    ## needed by functions oa.maxGR and oa.maxGRmin34
    IDname <- gsub("\"","",deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)){ 
    if (!exists(IDname)) 
          ID <- eval(parse(text=paste("oa.design(",IDname,")")))
    else if (is.character(ID)) 
          ID <- eval(parse(text=paste("oa.design(",IDname,")")))
    }
    GR <- 3
    hilf <- P3.3(ID, digits=digits^2, rela=TRUE)
    if (max(hilf[,1])==0){
       GR <- 4
       if (ncol(ID)>=4){
       hilf <- P4.4(ID, digits=digits^2, rela=TRUE)
          if (max(hilf[,1])==0){
             GR <- 5
          }}
       }
    RPFT <- NULL
    if (GR < ncol(ID) + 1){ 
      RPFT <- hilf
      GR <- round(GR + 1-sqrt(max(hilf[,1])),digits)
      if (GR==6) GR <- 5 ## maximum set to 5, as more has not been checked
      }
    else GR <- Inf
    list(GR=GR, RPFT=RPFT)
}