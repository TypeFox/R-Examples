read.convert_1_2 <- function(proc.length)
{

version<-read.table("wave_cor_mat_level_1.txt",skip=0,nrows=1)

tmp<-read.table("wave_cor_mat_level_1.txt",skip=1,nrows=1)

if(tmp[[2]]!="Correlation") stop("The files are not of class correlation matrix")

tmp<-read.table("wave_cor_mat_level_1.txt",skip=2,nrows=1)
method<-tmp[[1]]
wf<-tmp[[2]]
boundary<-tmp[[3]]
n.levels<-tmp[[4]]

wave.cor.list <- vector("list", (3*n.levels))
    names(wave.cor.list) <- c(paste("d", 1:n.levels, sep = ""), paste("lowerd", 1:n.levels, sep = ""),paste("upperd", 1:n.levels, sep = ""))

class(wave.cor.list)<-"Wave Correlation"
attr(wave.cor.list, "version") <- 2
attr(wave.cor.list, "method") <- method
attr(wave.cor.list, "wavelet") <- wf
attr(wave.cor.list, "boundary") <- boundary
attr(wave.cor.list, "proc.length") <- proc.length


for(i in 1:n.levels){

name.txt<-paste("wave_cor_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")


wave.cor.list[[i]]<-as.matrix(read.table(name.txt,skip=3))

name.txt<-paste("wave_cor_lower_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

wave.cor.list[[i+4]]<-as.matrix(read.table(name.txt,skip=3))

name.txt<-paste("wave_cor_upper_mat_level_",i,sep="")
name.txt<-paste(name.txt,".txt",sep="")

wave.cor.list[[i+8]]<-as.matrix(read.table(name.txt,skip=3))

}

return(wave.cor.list)

}

