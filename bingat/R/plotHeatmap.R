plotHeatmap <-
function(data, type, names, ...){
if(missing(data) || missing(type))
stop("data and/or type is missing.")

#Take only the first column of data if it is multi columned
if(class(data) == "data.frame" || class(data) == "matrix")
data <- data[,1]

mat <- vec2mat(data, type)

if(missing(names))
names <- 1:ncol(mat)
colnames(mat) <- names
rownames(mat) <- names

colfunc <- colorRampPalette(c("blue", "grey"))

gplots::heatmap.2(mat, symm=TRUE, Rowv=NA, dendrogram="none", trace="none", col=colfunc(10), ...)
}
