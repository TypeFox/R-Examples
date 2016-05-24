colormap<-function(design,mod) {
##################### Inputs ###########################################
# des - a data frame containing columns of the numeric factor levels
# mod - the model for the color plot of correlations
#    1 = Linear model containing only the terms in the data frame
#    2 = Linear model plus two factor interactions
#    3 = Linear model plus 2 and 3 factor interactions
#    4 = Linear model plus 2, 3 and 4 factor interactions
########################################################################
y<-runif(nrow(design),0,1)
if(mod==1) {test <- model.matrix(lm(y~(.),data=design))}
if(mod==2) {test <- model.matrix(lm(y~(.)^2,data=design))}
if(mod==3) {test <- model.matrix(lm(y~(.)^3,data=design))}
if(mod==4) {test <- model.matrix(lm(y~(.)^4,data=design))}
names<-colnames(test)
names<-gsub(':','',names)
names<-gsub('1','',names)
colnames(test)<-names
cmas<-cor(test[,ncol(test):2])
cmas<-cmas[c((ncol(cmas)):1), ]
rgb.palette <- colorRampPalette(c("white", "black"), space = "rgb")
levelplot(cmas, main="Color map of correlations", xlab="", ylab="", col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01),scales=list(x=list(rot=90)))
                          }
