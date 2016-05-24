######################
# MetaArray S4 class #
######################

myRep <- representation(GEDM = "list",
 				clinical = "list",
				datanames = "character")
 setClass("MetaArray", representation = myRep)

 setGeneric("print")
setMethod("print","MetaArray", function(x, ...) {
for (i in 1:length(x@GEDM))
{
cat("Dataset", x@datanames[i], " containing" , dim(x@GEDM[[i]])[1], "probes and ", dim(x@GEDM[[i]])[2], "samples. \n" )
cat("Sumarization of samples: \n")
cl<-x@clinical[[i]]
for (j in 1:ncol(cl)) {
cat(colnames(cl)[j],"\n")
cat( paste(rep("-",nchar(colnames(cl)[j])), collapse=""),"\n")
print(summary(cl[,j])) 
}
cat("\n")
}
})

setMethod("show","MetaArray", function(object) {
for (i in 1:length(object@GEDM))
{
cat("Dataset", object@datanames[i], " containing" , dim(object@GEDM[[i]])[1], "probes and ", dim(object@GEDM[[i]])[2], "samples. \n" )
cat("Sumarization of samples: \n")
cl<-object@clinical[[i]]
for (j in 1:ncol(cl)) {
cat(colnames(cl)[j],"\n")
cat( paste(rep("-",nchar(colnames(cl)[j])), collapse=""),"\n")
print(summary(cl[,j])) 
}
cat("\n")
}
})

 setGeneric("as.list")
setMethod("as.list", "MetaArray", function(x){
N<-length(GEDM(x))
l<-list()
for (i in 1:N) {
l[[i]]<-list(GEDM=GEDM(x)[[i]], clinical=clinical(x)[[i]], dataname=datanames(x)[[i]])
}
return(l)
}
)


clinical <- function(object) {return(object@clinical)} 
GEDM <- function(object) {return(object@GEDM)}      
datanames <-function(object) {return(object@datanames)}  

 setGeneric("GEDM<-", function(object, value) standardGeneric("GEDM<-"))
 setReplaceMethod("GEDM", "MetaArray", function(object, value) {
 object@GEDM <- value
 object
 })
 setGeneric("clinical<-", function(object, value) standardGeneric("clinical<-"))
 setReplaceMethod("clinical", "MetaArray", function(object, value) {
 object@clinical <- value
 object
 })

 setGeneric("datanames<-", function(object, value) standardGeneric("datanames<-"))
 setReplaceMethod("datanames", "MetaArray", function(object, value) {
 object@datanames <- value
 object
 })




setMethod("plot","MetaArray", function(x, col= c("red", "blue", "green"),...) {
	CL<-clinical(x)
datname<-datanames(x)

N<-length(CL)
TT<-clinical.sum(x)$absolute

totalSamples<-sum(sapply(CL,nrow))
param<-names(TT)
Nparam<-length(param)
Nchar=max(nchar(datanames(x)))

grid.newpage()
vplay<-grid.layout(N+2,Nparam+2, 
	widths =unit(c(Nchar-3, 1, rep(1,Nparam)), c("char", "char", rep("null", Nparam))))

vp<-viewport(layout=vplay)

pushViewport(vp)
#Names and Sample Sizes
for (i in 1:N)
{
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=i, margins=c(0,0,0,0)))
grid.text(datanames(x)[i])
popViewport(1)
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=i, margins=c(0,0,0,0)))
grid.text(nrow(clinical(x)[[i]]))
popViewport(1)
}

#Clinical Parameters
con<-unlist(lapply(TT, function(x) "Mean" %in% colnames(x)))
for (j in 1:Nparam){

if (!con[j])
{

  for (i in 1:N)
  {
  pushViewport(viewport(layout.pos.col=2+j, layout.pos.row=i))
  if (param[j] %in% colnames(CL[[i]]) ){
    X<-prop.table(t(TT[[param[j]]][i,]))
    for (k in 1:length(X)) 
    grid.rect(x = unit(0.05+0.9*sum(X[1,1:k]), "npc"), y = unit(0.05, "npc"),
              width = unit(0.9*X[1,k], "npc"), height = unit(0.9, "npc"),
              just = c("right", "bottom"), gp=gpar(fill=col[k]))
    } 
  popViewport(1)
  }

} else {
for (i in 1:(N))
  {
   if (param[j] %in% colnames(CL[[i]])){
    b<-TT[[param[j]]][i,c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")]
    bb<-TT[[param[j]]][N+1,c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")]
    pushViewport(viewport(layout.pos.col=2+j, layout.pos.row=i,xscale=c(bb[1]-10,bb[5]+10)))
     grid.lines(x=unit(c(b[1], b[5]),"native"), y=unit(c(0.5, 0.5),"npc"))
     grid.rect(x=unit(b[2], "native"), width=unit(b[4]-b[2],"native"), height=unit(0.6, "npc"), 
      just=c("left", "center"), gp=gpar(fill=col))
     grid.rect(x=unit(b[3],"native"), width=unit(1, "mm"), height=unit(0.6,"npc"),
      just=c("right","center"),gp=gpar(fill="black"))

     grid.lines(x=unit(b[1],"native"), y=unit(c(0.4,0.6),"npc") )
     grid.lines(x=unit(b[5],"native"), y=unit(c(0.4,0.6),"npc") )
    popViewport(1)
    } 
  }
}

}

#Overall 
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=N+1, margins=c(0,0,0,0)))
grid.text("Total")
popViewport(1)
pushViewport(plotViewport(layout.pos.col=2, layout.pos.row=N+1, margins=c(0,0,0,0)))
grid.text(sum(sapply(CL, nrow)))
popViewport(1)

for (j in 1:Nparam){

if (!con[j])
{

  pushViewport(viewport(layout.pos.col=2+j, layout.pos.row=N+1))
    X<-prop.table(t(TT[[param[j]]][N+1,]))
    for (k in 1:length(X)) 
    grid.rect(x = unit(0.05+0.9*sum(X[1,1:k]), "npc"), y = unit(0.05, "npc"),
              width = unit(0.9*X[1,k], "npc"), height = unit(0.9, "npc"),
              just = c("right", "bottom"), gp=gpar(fill=col[k]))
  
  popViewport(1)
  
} else {
    bb<-TT[[param[j]]][N,c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")]
    pushViewport(viewport(layout.pos.col=2+j, layout.pos.row=N+1,xscale=c(bb[1]-10,bb[5]+10)))
     grid.lines(x=unit(c(bb[1], bb[5]),"native"), y=unit(c(0.5, 0.5),"npc"))
     grid.rect(x=unit(bb[2], "native"), width=unit(bb[4]-bb[2],"native"), height=unit(0.6, "npc"), 
      just=c("left", "center"), gp=gpar(fill=col))
     grid.rect(x=unit(bb[3],"native"), width=unit(1, "mm"), height=unit(0.6,"npc"),
      just=c("right","center"),gp=gpar(fill="black"))

     grid.lines(x=unit(bb[1],"native"), y=unit(c(0.4,0.6),"npc") )
     grid.lines(x=unit(bb[5],"native"), y=unit(c(0.4,0.6),"npc") )
     grid.xaxis()
    popViewport(1)
    
}

}

#Legend
is.odd <- function(x) x %% 2 != 0 
for (j in 1:Nparam){
pushViewport(viewport(layout.pos.col=2+j, layout.pos.row=N+2, layout=vplay))
labels=colnames(TT[[param[j]]])
nkeys=length(labels)
if (con[j]) grid.text(param[j], x = unit(0.5, "npc"), y = unit(0.1, "npc"),just=c("centre","bottom")) else {
for (k in 1:nkeys)
{
if (is.odd(k)) {
grid.rect(x = unit(0.1, "npc"), y = unit(1-0.3*ceiling(k/2) , "npc"),
               width = unit(0.1, "npc"), height = unit(0.1, "npc"),
               just = c("centre", "center"), gp=gpar(fill=col[k])) 
grid.text(labels[k], x = unit(0.2, "npc"), y = unit(1-0.3*ceiling(k/2), "npc"),just=c("left","center"))
} else {
grid.rect(x = unit(0.6, "npc"), y = unit(1-0.3*(ceiling(k/2)), "npc"),
               width = unit(0.1, "npc"), height = unit(0.1, "npc"),
               just = c("centre", "center"), gp=gpar(fill=col[k])) 
grid.text(labels[k], x = unit(0.7, "npc"), y = unit(1-0.3*(ceiling(k/2)), "npc"),just=c("left","center"))
}
}
grid.text(param[j], x = unit(0.5, "npc"), y = unit(0.1, "npc"),just=c("centre","bottom"))
}
popViewport(1)
}


}
)
