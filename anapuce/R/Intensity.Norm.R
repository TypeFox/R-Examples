Intensity.Norm <-
function (fileIN = "resNorm.txt", n = 3, ind.array = NULL, name.A = "A", 
    name.M = "M.norm", sep="\t", center=FALSE,log.transf=TRUE, ...) 
{
    
    if (is.character(fileIN)) {
    res <- read.table(fileIN, header = TRUE, sep = sep, ...)
    } else {
        res <- fileIN
    } 
     all.array<- sub(name.M,"",names(res)[grep(name.M,names(res))])
    if (is.null(ind.array))
       ind.array<-all.array
    else {
       if(length(setdiff(ind.array,all.array))!=0)
         cat(paste(name.M,setdiff(ind.array,all.array),sep="")," not found and not considered \n") 
       ind.array<-intersect(ind.array,all.array)   
        }
    indA <- which(names(res) %in% paste(name.A,ind.array,sep=""))
    indM <- which(names(res) %in% paste(name.M,ind.array,sep=""))
    
    if (center){
        ## Modif parallelisation
        ##MeanArray=apply(as.matrix(res[, indA]),2,mean,na.rm=TRUE)
        MeanArray <- colMeans(as.matrix(res[, indA]),na.rm=TRUE)
        MeanAll = mean(MeanArray)
        A = apply(as.matrix(1:length(indA)),1,FUN=function(x) res[, indA[x]]-MeanArray[x]+MeanAll)
        
    } else {
        A = res[,indA]
    }  
     ind=rbind(indA,indM)
     Red = apply(ind,2,FUN=function(x) (res[, x[1]] + (0.5 * res[, x[2]])) )
     Green = apply(ind,2,FUN=function(x) (res[, x[1]] - (0.5 * res[, x[2]])) )
     if(log.transf)
      Int = data.frame(res[,1:n],Red,Green)
     else 
      Int = data.frame(res[,1:n],2^(Red),2^(Green))
     names(Int)=c(names(res)[1:n],paste("Red.norm",ind.array,sep=""),paste("Green.norm",ind.array,sep=""))
   
      Int = Int[,c(1:n,match( paste(c("Red.norm","Green.norm"),rep(ind.array,each=2),sep=""),names(Int)))]
 
    invisible(Int)
    # (c) 2007 Institut National de la Recherche Agronomique

    
}

