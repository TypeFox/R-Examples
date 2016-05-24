	#make block design matrix 
makeBlkDesMatrix <-
function(design.df, blkTerm){
isFactorNameNumeric <-
function(levels) !as.logical( length(grep("[A-Z]|[a-z]", levels)) )

makeDesignMatrix <-
function(nRows, design.df, col){    
if(grepl(":", col)){ 
factor <-as.factor( apply(design.df[,unlist(strsplit(col, ":"))], 1, function(x) paste(x, collapse =".")))
}else{
factor <- as.factor(design.df[,col])
}  

facName <- col
nCols <- nlevels(factor)

Z <- matrix(0, nrow=nRows, ncol=nCols)
Z[cbind(1:nRows, match(c(factor), 1:nCols))] <- 1
if(isFactorNameNumeric(levels(factor))){ 
colNames <- paste(facName, 1:nCols,sep="")
}else{                                    
colNames <- levels(factor)
}

dimnames(Z) <- list(1:nRows, colNames)
return(Z)
}

# design.df = data.frame containing design
# blkTerm = block terms

n <- length(blkTerm)
nRows <- nrow(design.df)
Z <- list(NULL)
Z[[1]] <- diag(nrow(design.df))

for(i in 2:(n+1)){   
Z[[i]] <- makeDesignMatrix(nRows=nRows, design.df=design.df, col=blkTerm[i-1]) 
}   

names(Z) <- c("e", blkTerm)  
return(Z)
}
