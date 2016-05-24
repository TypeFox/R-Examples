## ----, message=FALSE-----------------------------------------------------
library(gapmap)

## ------------------------------------------------------------------------
set.seed(1234)
x <- rnorm(10, mean=rep(1:5, each=2), sd=0.4)
y <- rnorm(10, mean=rep(c(1,2), each=5), sd=0.4)
dataFrame <- data.frame(x=x, y=y, row.names=c(1:10))
#calculate distance matrix. default is Euclidean distance
distxy <- dist(dataFrame)
#perform hierarchical clustering. default is complete linkage.
hc <- hclust(distxy)
dend <- as.dendrogram(hc)

## ----, fig.width= 6.5, fig.height=6--------------------------------------
grey_scale =c("#333333", "#5C5C5C", "#757575", "#8A8A8A", "#9B9B9B", "#AAAAAA", "#B8B8B8", "#C5C5C5", "#D0D0D0", "#DBDBDB", "#E6E6E6")
gapmap(m = as.matrix(distxy), d_row= rev(dend), d_col=dend, col = grey_scale)

## ----, fig.width= 6.5, fig.height=6--------------------------------------
gapmap(m = as.matrix(distxy), d_row= rev(dend), d_col=dend,  mode = "quantitative", mapping="linear", col = grey_scale)

## ----, echo=FALSE, fig.width= 6.5, fig.height=6--------------------------
distances <-seq(0, 5, 0.1)
data <- data.frame(distance=distances)
s <- 0.5
l <- data
e <- data
for(i in 1:nrow(data)){
  dist <- data$distance[i]
  linear <- map(dist, 0, 5, 0, 1)
  exp <-  map.exp(dist, 0, 5, 0, 1, scale = s)
  #print(paste0("dist =", dist," linear=",linear, " exp=", exp))
  l$gap[i] = linear
  e$gap[i] = exp
  l$type[i] = "linear"
  e$type[i] = "exponential"
}
gaps <- rbind(l, e)
ggplot(gaps, aes(x=gap, y=distance, group=type)) + geom_line(aes(color=type))+theme_bw()+ theme(legend.position= c(0.9,0.1))

## ----, echo=FALSE, fig.width= 6.5, fig.height=6--------------------------
scales <- seq(0.1, 3, 0.3)
distances <-seq(0, 5, 0.1)
D = data.frame()
for(j in 1:length(scales)){
  s  <- scales[j]
  data <- data.frame(distance=distances)
  for(i in 1:nrow(data)){
    dist <- data$distance[i]
    exp <-  map.exp(dist, 0, 5, 0, 1, scale = s)
    data$gap[i] = exp
    data$scale[i] = s
  }
  D <- rbind(D, data)
}
labels = data.frame()
for(j in 1:length(scales)){
  a = 0
  b = 5
  c = 0
  d = 1
  y = 0.4 # x position  
  x = scales[j]
  v = a + ((y/(d-c))^x) *(b-a) 
  labels <- rbind(labels, data.frame(scale=x, distance =v, gap=y))
}
ggplot() + geom_line(data=D, aes(x=gap, y=distance, group=scale), color="#56B1F7")+ scale_y_continuous(limits = c(0,5))+
  geom_text(data= labels, aes(x=gap,y=distance, label=scale), hjust=-0.2, vjust=0) +
  geom_point(data= labels, aes(x=gap,y=distance)) +
  theme_bw() + theme(legend.position="none") 

## ----, fig.width= 6.5, fig.height=6--------------------------------------
gapmap(m = as.matrix(distxy), d_row= rev(dend), d_col=dend,  mode = "threshold", row_threshold = 2, col_threshold = 2, col = grey_scale)

## ----, fig.width= 6.5, fig.height=6--------------------------------------
library(dendsort);
gapmap(m = as.matrix(distxy), d_row= rev(dendsort(dend)), d_col=dendsort(dend),  mode = "quantitative", col = grey_scale)

