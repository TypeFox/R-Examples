`vit` <-
function(id="", occasion="", score="", Data=NULL, group=NULL, subset.ids=NULL, pct.rand=NULL, 
number.rand = NULL, All.in.One=TRUE, ylab=NULL, xlab=NULL, same.scales=TRUE, plot.points=TRUE,
save.pdf=FALSE, save.eps=FALSE, save.jpg=FALSE, file="", layout = c(3,3), col=NULL, pch=16, 
cex=.7, ...)

{

if(is.null(Data)) stop("You need to specify the data set with the \'Data\' argument.")

if(id=="") stop("You need to specify ID variable with the \'id\' argument.")
if(score=="") stop("You need to specify score variable with the \'score\' argument.")
if(occasion=="") stop("You need to specify time variable with the \'occasion\' argument.")

if(file=="") 
{
file <- "vit"
no.file.name <- TRUE
}
if(save.pdf==TRUE) 
{
if(save.eps==TRUE) stop("Specify one file format at one time")
if(save.jpg==TRUE) stop("Specify one file format at one time")
pdf(file = paste(file,".pdf",sep=""))
}

if(save.eps==TRUE) 
{
if(save.jpg==TRUE) stop("Specify one file format at one time")
postscript(file = paste(file,".eps",sep=""))
}

if(save.jpg==TRUE) jpeg(filename = paste(file,"%d.jpg",sep=""),width = 640, height = 550)

# '%d' is used for plotting more than one page on one of these devices, 
# it retains files with the sequence numbers.

if(is.null(ylab)) ylab=score
if(is.null(xlab)) xlab=occasion
if(!is.null(group)) group.title=group

if(is.data.frame(Data))
{
id <- Data[,which(names(Data)==id)]
occasion <- Data[,which(names(Data)==occasion)]
score <- Data[,which(names(Data)==score)]
if(!is.null(group)) 
{
group <- Data[,which(names(Data)==group)]
G <- length(unique(group))
}
}

if(!is.data.frame(Data))
{
id <- Data[,(colnames(Data)==id)]
occasion <- Data[,(colnames(Data)==occasion)]
score <- Data[,(colnames(Data)==score)]
if(!is.null(group)) 
{
group <- Data[,(colnames(Data)==group)]
G <- length(unique(group))
}
}

if(!is.null(subset.ids))
{

if(!is.null(pct.rand)) stop("Since \'subset.ids\' was specified, do not also specify \'pct.random\'.")
if(!is.null(number.rand)) stop("Since \'subset.ids\' was specified, do not also specify \'number.rand\'.")

Rows.to.Select <- matrix(NA, length(id), length(subset.ids))
for(i in 1:length(subset.ids))
{
Rows.to.Select[,i] <- id==subset.ids[i]
}
Rows.to.Select <- apply(Rows.to.Select, MARGIN=1, sum)==1

id <- id[Rows.to.Select]  
occasion <- occasion[Rows.to.Select] 
score <- score[Rows.to.Select]
if(!is.null(group)) 
{
group <- group[Rows.to.Select]
G <- length(unique(group))
}
}

if(!is.null(number.rand))
{
if(!is.null(pct.rand)) stop("Since \'number.rand\' was specified, do not also specify \'pct.rand\'.")
}


if(is.null(subset.ids))
{
if(is.null(number.rand)) 
{
if(is.null(pct.rand)) pct.rand <- 1
if(pct.rand <=1 & pct.rand>0) number.rand <- ceiling(length(unique(id))*pct.rand)
if(pct.rand >1 & pct.rand<=100) number.rand <- ceiling(length(unique(id))*pct.rand/100)
}

rand.samp <- sample(unique(id), number.rand, replace=FALSE)

Data.From.Rand <- matrix(NA, length(id), length(rand.samp))
for(i in 1:length(rand.samp))
{
Data.From.Rand[,i] <- id==rand.samp[i]
}

Rows.to.Select <- apply(Data.From.Rand, MARGIN=1, sum)==1
id <- id[Rows.to.Select]  
occasion <- occasion[Rows.to.Select] 
score <- score[Rows.to.Select]
if(!is.null(group)) group <- group[Rows.to.Select]
}

if(!is.null(group)) 
{
group <- unique(cbind(as.character(id), as.character(group)))[,2]
selected.G <- length(unique(group))
}

ID <- unique(id)
N <- length(ID)


# Up until here, the function has only organized that data that is to be used. 
####################################################################################

# Optionally sets up the limits of the axes automatically.

ylim <- NULL
xlim <- NULL

if(same.scales==TRUE) 
{
ylim <- c(min(score),max(score))
xlim <- c(min(occasion),max(occasion))
}

####################################################################################

# Draws the plotting region.
if(All.in.One==FALSE)
{
par(mfrow=layout, mar = c(4, 3, 2, 1), mgp = c(2, .5, 0), ...)
####################################################################################
# The following optionally plots the values of the observed data.
if(N > layout[1]*layout[2])
{
if(save.pdf==FALSE & save.eps==FALSE & save.jpg==FALSE) 
{
print(paste("In order to display multiple graphic windows (e.g., if there are more than", layout[1]*layout[2],"participants), click \'Recording\' in the \'History\' menu of an open plot window. Using \'Page Up\' and \'Page Down\' will scroll through the graphic windows."))
}
}

if(!is.null(group)) 
{
color <- matrix(NA, N, 1)
if(G!=length(col)&&!is.null(col))
{
print("number of groups does not equal to number of colors specified, only first color was used") 
color <- rep(col[1],N)
}
if(G==length(col))
{
for (j in 1:selected.G)
{
for (i in 1:N)
{
if(group[i]==unique(group)[j]) color[i] <- col[j]
}
}
}
if(is.null(col))
{
for (j in 1:selected.G)
{
for (i in 1:N)
{
if(group[i]==unique(group)[j]) color[i] <- j
}
}
}

for(i in 1:N)
{

plot(occasion[id==ID[i]], score[id==ID[i]], ylim=ylim, xlim=xlim, ylab=ylab, xlab=xlab, type="n", 
     main=paste("(ID = ",ID[i],", " ,group.title," = ",group[i],")",sep=""),
     
     
     #substitute({" (ID" == ID.i}* {", Group" == group.i} * {")"}, 
     #list(ID.i = ID[i],group.i = group[i])),
     cex.main=.9, cex.lab=.9, ...)


if(plot.points==TRUE) points(occasion[id==ID[i]], score[id==ID[i]], col=color[i],pch=pch,cex=cex, ...)
lines(occasion[id==ID[i]], score[id==ID[i]], col=color[i],  ...)
}
}

if(is.null(group)) 
{
if(is.null(col)) col <- 1
for(i in 1:N)
{

plot(occasion[id==ID[i]], score[id==ID[i]], ylim=ylim, xlim=xlim, ylab=ylab, xlab=xlab, type="n", 
     main=paste("(ID = ",ID[i],")",sep=""),

     #substitute({" (ID" == ID.i} * {")"}, list(ID.i = as.character(ID)[i])),
     cex.main=.9, cex.lab=.9, ...)

if(plot.points==TRUE) points(occasion[id==ID[i]], score[id==ID[i]],col=col[1],pch=pch,cex=cex, ...)
lines(occasion[id==ID[i]], score[id==ID[i]],col=col[1], ...)
}
}
}

####################################################################################
if(All.in.One==TRUE)
{

if(!is.null(group)) 
{
color <- matrix(NA, N, 1)
if(G!=length(col)&&!is.null(col))
{
print("number of groups does not equal to number of colors specified, only first color was used") 
color <- rep(col[1],N)
}
if(G==length(col))
{
for (j in 1:selected.G)
{
for (i in 1:N)
{
if(group[i]==unique(group)[j]) color[i] <- col[j]
}
}
}
if(is.null(col))
{
for (j in 1:selected.G)
{
for (i in 1:N)
{
if(group[i]==unique(group)[j]) color[i] <- j
}
}
}

plot(0, ylim=c(min(score),max(score)), xlim=c(min(occasion),max(occasion)), ylab=ylab, xlab=xlab, font.lab=3, type="n", ...)

for(i in 1:N)

{
if(plot.points==TRUE) points(occasion[id==ID[i]], score[id==ID[i]], col=color[i],pch=pch,cex=cex, ...)
lines(occasion[id==ID[i]], score[id==ID[i]], col=color[i],  ...)
}
}


if(is.null(group)) 
{
plot(0, ylim=c(min(score),max(score)), xlim=c(min(occasion),max(occasion)), ylab=ylab, xlab=xlab, font.lab=3, type="n", ...)
if(is.null(col)) col <- 1
for(i in 1:N)

{
if(plot.points==TRUE) points(occasion[id==ID[i]], score[id==ID[i]],col=col[1],pch=pch,cex=cex, ...)
lines(occasion[id==ID[i]], score[id==ID[i]],col=col[1], ...)
}
}


}

####################################################################################
par(mfrow=c(1,1))

if(save.pdf==TRUE)
{
dev.off()
if(no.file.name==TRUE) print(paste("\'vit.pdf\' file saved at the directory",getwd()))
}

if(save.eps==TRUE)
{
dev.off()
if(no.file.name==TRUE) print(paste("\'vit.eps\' file saved at the directory",getwd()))
}

if(save.jpg==TRUE)
{
dev.off()
if(no.file.name==TRUE) print(paste("\'vit.jpg\' file(s) saved at the directory",getwd()))
}

}

