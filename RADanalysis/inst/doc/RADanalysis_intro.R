## ------------------------------------------------------------------------
#install.packages("RADanalysis")
library(RADanalysis)

## ---- echo=TRUE, results='asis',fig.cap="rank"---------------------------
rad <- data.frame(rank = 1:100, abundance = round(1./(1:100) * 100))
knitr::kable(head(rad, 10))

## ---- fig.show='hold',fig.width=5,fig.height=5,fig.cap="A typical Rank Abundance Distribution in linear-linear scale"----
plot(rad,xlab = "Rank",ylab = "Abundance",pch = 19,type = "b",lwd = 0.5)

## ---- fig.show='hold',fig.width=5,fig.height=5,fig.cap="A typical Rank Abundance Distribution in log-log scale"----
plot(rad,xlab = "Rank",ylab = "Abundance",pch = 19,log = "xy",type = "b",lwd = 0.5)

## ---- echo=TRUE, results='asis',fig.cap="rank"---------------------------
library(RADanalysis)
data("gut_otu_table")

colnames(gut_otu_table) <- c(paste("OTU_",seq(1:ncol(gut_otu_table)),sep = ""))
knitr::kable(gut_otu_table[,1:10])

## ---- fig.show='hold',fig.width=6,fig.height=6,fig.cap="Rank Abundance Distributions of gut microbiome data",warning=FALSE----
data("gut_otu_table")
rads <- gut_otu_table

#plot original rads
line_cols <- c("green","red","blue")
#to specify different stages of subjects
sample_classes <- c(1,1,1,1,2,2,3,3,1,1,2,3,3,1,1,2,3,3)
plot(1,xlim = c(1,2000),ylim = c(1,20000),col = "white",log  = "xy",axes = F,
     xlab = "Rank",ylab = "Abundance",main = "")
sfsmisc::eaxis(side = 1,at = c(1,10,100,1000))
sfsmisc::eaxis(side = 2)
for(i in 1:nrow(rads)){
    temp <- sort(rads[i,],decreasing = TRUE)
    temp <- temp[temp>0]
    lines(x = temp,lwd = 2,col = line_cols[sample_classes[i]])
}
legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),
       col = line_cols,lwd = 3)

## ---- fig.show='hold',fig.width=8,fig.height=8,fig.cap="Normalizing RADs using MaxRank normalization method with RADnormalization()",warning=FALSE----
data(gut_otu_table)
rads <- gut_otu_table

original_rad <- sort(rads[1,],decreasing = TRUE)
#removing zeros
original_rad <- original_rad[original_rad > 0]
plot(original_rad,ylim = c(1,max(original_rad)),log = "xy", xlab = "Rank",
     ylab = "Abundance", main = "",pch = 19,type = "b",cex = 0.5)

norm_rad <- RADnormalization(input = rads[1,],max_rank = 400,average_over = 50)
points(x = norm_rad$norm_rad * sum(norm_rad$norm_rad_count[1,]) ,pch = 19,cex = 1, type = "l",
       col = "blue",lwd = 4)
points(x = norm_rad$norm_rad_count[1,],pch = 19,cex = 1, type = "l",col = "red",lwd = 3)
points(x = norm_rad$norm_rad_count[2,],pch = 19,cex = 1, type = "l",col = "green",lwd = 3)
legend("bottomleft",legend = c("original RAD","possible norm rad","possible norm rad",
                    paste("nrad averaged over 50 realizations, times", 
                          sum(norm_rad$norm_rad_count[1,]))),
                    col = c("black","red","green","blue"),lwd = 2,bty = "n")

## ---- fig.show='hold',fig.width=7,fig.height=3,fig.cap="Richness distribution in gut_otu_table",warning=FALSE----

data("gut_otu_table")
rads <- gut_otu_table

richness <- sapply(X = 1:nrow(rads), function(i) length(which(rads[i,] > 0)))
boxplot(richness,horizontal = T,xlab = "Richness")
quantile(richness)

## ---- fig.show='hold',fig.width=7,fig.height=7,fig.cap="NRADs of gut_otu_table with $R=400$ and averaged over 10 possible NRADs",warning=FALSE----
data("gut_otu_table")
rads <- gut_otu_table

#plot original rads
line_cols <- c("green","red","blue")
sample_classes <- c(1,1,1,1,2,2,3,3,1,1,2,3,3,1,1,2,3,3)

#Normalization
nrads <- RADnormalization_matrix(input = rads,max_rank = 400,average_over = 10,
                                 sample_in_row = TRUE,verbose = FALSE)
nrads <- nrads$norm_matrix

plot(1,xlim = c(1,400),ylim = c(4e-5,1),col = "white",log  = "xy",
     xlab = "Rank",ylab = "Abundance",
     main = "")
for(i in 1:nrow(nrads)){
    lines(x = nrads[i,],lwd = 2,col = line_cols[sample_classes[i]])
}
legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),
       col = line_cols,lwd = 3)

## ---- fig.show='hold',fig.width=7,fig.height=7,warning=FALSE,fig.cap="Representative NRADs of three states of samples"----
data("gut_nrads")
nrads <- gut_nrads
nrads <- nrads$norm_matrix

line_cols <- c("green","red","blue")
sample_classes <- c(1,1,1,1,2,2,3,3,1,1,2,3,3,1,1,2,3,3)
maxrank <- 400


#plot nrads
plot(1e10,xlim = c(1,maxrank),ylim = c(2e-5,1),log="xy",
     xlab = "rank",ylab = "abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,10,100,1000,10000))
sfsmisc::eaxis(side = 2,at = c(1e-4,1e-3,1e-2,1e-1,1),las = 0)

for(i in 1:nrow(nrads)){
  points(nrads[i,],type = 'l',col = line_cols[sample_classes[i]],lwd = 0.8)
}
#plot confidence intervals of representative nrads 
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 1),
                      plot = TRUE,confidence = 0.9,with_conf = TRUE,
                      col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 2),
                      plot = TRUE,confidence = 0.9,with_conf = TRUE,
                      col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 3),
                      plot = TRUE,confidence = 0.9,with_conf = TRUE,
                      col = scales::alpha(line_cols[3],0.5),border = NA)
#plot representative nrads
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 1),
                      plot = TRUE,with_conf = FALSE,
                      col = scales::alpha(line_cols[1],0.8),lwd = 4)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 2),
                      plot = TRUE,with_conf = FALSE,
                      col = scales::alpha(line_cols[2],0.8),lwd = 4)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == 3),
                      plot = TRUE,with_conf = FALSE,
                      col = scales::alpha(line_cols[3],0.8),lwd = 4)
legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),
       col = line_cols,lwd = 3)

## ---- fig.show='hold',fig.width=7,fig.height=7,warning=FALSE,fig.cap="Multi dimensional scaling of gut data using manhattan distance. Large points and errorbars are representative points and standard error of the mean for each group."----
data("gut_nrads")
nrads <- gut_nrads
nrads <- nrads$norm_matrix

line_cols <- c("green","red","blue")
sample_classes <- c(1,1,1,1,2,2,3,3,1,1,2,3,3,1,1,2,3,3)
maxrank <- 400

#distance matrix using manhattan distance
d <- dist(x = nrads,method = "manhattan")
#ordination using classical multi-dimensional scaling
mds <- cmdscale(d = d,k = 5,eig = TRUE)

#plot the points 
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS plot with representative points \n of each group and error bars")

#add the representative points wit erorr bar to the previous plot
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)

legend("bottomleft",bty = "n",legend = c("pre Cp","under Cp","post Cp"),
       col = line_cols,pch = 19)

