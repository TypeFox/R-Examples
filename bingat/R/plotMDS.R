plotMDS <-
function(data, groups, estGstar=TRUE, paired=FALSE, returnCoords=FALSE, ...){
if(missing(data))
stop("data is missing.")

if(missing(groups))
groups <- rep(1, ncol(data))

if(length(groups) != ncol(data))
stop("groups must be the same length as the data set.")

if(min(groups) <= 0)
groups <- groups - min(groups) + 1

dlen <- ncol(data)
uni <- unique(groups)
mycol <- groups + 1
mycex <- rep(1, length(groups))
numcols <- max(length(uni), 2, max(mycol))

if(paired){
if(length(uni) != 2 || sum(groups==uni[1]) != sum(groups==uni[2]))
stop("There must be exactly TWO groups of the same length when using paired==TRUE.")
groupSize <- ncol(data)/2
}

if(estGstar){
gss <- NULL
if(length(uni) > 1){
for(i in 1:length(uni)){
if(sum(groups==uni[i]) == 1)
next
gss <- cbind(gss, estGStar(data[,groups==uni[i]]))
colnames(gss)[ncol(gss)] <- paste("Group", i, "Gstar")
mycol <- c(mycol, max(mycol)+1)
}
}
gss <- cbind(gss, estGStar(data))
numcols <- numcols + ncol(gss)
colnames(gss)[ncol(gss)] <- paste("Combined Gstar")
mycol <- c(mycol, 1)
mycex <- c(mycex, rep(1.5, ncol(gss)))
data <- cbind(data, gss)
}
palette(c("black", rainbow(numcols)))

loc <- cmdscale(dist(t(data), method="binary"), k=2) 
plot(loc, xlab="MDS 1", ylab="MDS 2", pch=16, col=mycol, cex=mycex, ...)

if(paired){
for(i in 1:groupSize)
segments(loc[i,1], loc[i,2], loc[i+groupSize,1], loc[i+groupSize,2], lty=3)

if(estGstar)
segments(loc[groupSize*2+1,1], loc[groupSize*2+1,2], loc[groupSize*2+2,1], loc[groupSize*2+2,2], lty=3) 
}

if(estGstar){
ptr <- nrow(loc) - ncol(gss)
for(i in 1:ncol(gss)){
text(loc[ptr+i,1], loc[ptr+i,2], rownames(loc)[ptr+i], pos=3)
}
}

if(returnCoords)
return(loc)
}
