polys.setup <- function(object){

#require("sp")

#requireNamespace("sp", quietly = TRUE)


obj <- fortify(object)
OBJ <- obj[, c("long","lat","piece","id")] 

pr <- c((min(as.numeric(names(table(OBJ[,3])))) + 1) : (max(as.numeric(names(table(OBJ[,3]))))) )  

ind.pr <- sum(as.numeric(table(OBJ[,3])!=0))


if(ind.pr > 1){

lldo <- dim(OBJ)[1]

lldo <- lldo + sum(as.numeric(OBJ[2:lldo, 3] != OBJ[1:(lldo-1),3])) 


for(i in 1:lldo){ 
    
    if( OBJ[i,3] %in% pr && ( as.numeric(OBJ[i,3]) - as.numeric(OBJ[i-1,3]) ) == 1)  {

                 new <- as.numeric( c(NA, NA, OBJ[i,3], OBJ[i,4]) )
                 OBJ <- rbind(OBJ, new)[c( 1:(i-1) , dim(OBJ)[1]+1 , i:dim(OBJ)[1] ) , ]
                                                                                           }
}
    
    
}


OB  <- split(OBJ[,1:2], OBJ[,4]) 
llOB <- length(OB)

for(i in 1:llOB) OB[[i]]  <- data.matrix(OB[[i]]) 

list(polys = OB, names0 = object$NAME_0, names1 = object$NAME_1, names2 = object$NAME_2 )


}


