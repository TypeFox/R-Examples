getBranchSizes <-
function(tree, colors, divisions){
if(missing(tree) || class(tree) != "phylo")
stop("A valid tree is required.")

edgeColor <- NULL
edgeWidth <- NULL
edgeLength <- tree$edge.length

if(max(edgeLength) == 0){ #catch an all 0 tree
retData <- list(edgecol=rep(0, length(edgeLength)), edgewid=rep(0, length(edgeLength)))
return(retData)
}

if(missing(divisions))
divisions <- c(.1, 1, 10, 100, 1000, 10000)

divisions <- sort(divisions)

if(missing(colors))
colors <- c("red", "orange", "yellow", "green" , "cyan", "blue")

if(length(divisions) > (length(colors)+1)) #need more colors, dont care if more colors than divisons
colors <- c(colors, rep(colors[length(colors)], (length(divisions) - length(colors)-1)))

palette(colors)

for(i in 1:length(edgeLength)){ 
if(edgeLength[i] == 0){ #0 value so make it white
edgeColor <- append(edgeColor, 0)
edgeWidth <- append(edgeWidth, 0)
}else{
for(j in 1:length(divisions)){
if(edgeLength[i] <= divisions[j]){
cval <- j
len <- floor(4*(j+1)/length(divisions)) 
edgeColor <- append(edgeColor, cval)
edgeWidth <- append(edgeWidth, len)
break
}
}
}
}

retData <- list(edgecol=edgeColor, edgewid=edgeWidth)
return(retData)
}
