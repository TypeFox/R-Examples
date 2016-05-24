disc.mentr <-
function (data, varcon, out=c("symb","num")) 
{#Discretization using entropy with MDL as stopping criteria
    #Fixed by Edgar Acuna, November 2015
    n <- dim(data)[1]
    p <- dim(data)[2]
    data1 <- data
    for (j in varcon) {
        var <- varcon[j]
        sal <- discretevar(data, var, n, p)
        nparti <- sal[1]
        points <- sal[-1]
        if(out=="num"){
             if(nparti>1)
         {data1[, j] <- as.vector(cut(data[,j],nparti,labels=FALSE))}
 else{data1[,j]=1}}
else{
 if(nparti>1)
 {limites=c(-Inf,points,Inf)
     data1[, j] <- as.vector(cut(data[,j],limites))}
 else{data1[,j]="All"}}
}
    return(data1)
}
