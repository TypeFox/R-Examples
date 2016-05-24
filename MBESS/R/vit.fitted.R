`vit.fitted` <-
function(fit.Model, layout = c(3,3), ylab = "", xlab = "", pct.rand = NULL, number.rand = NULL,
                 subset.ids=NULL, same.scales = TRUE, save.pdf=FALSE, save.eps=FALSE, save.jpg=FALSE, file="", ...)
                 
{

###############################################################################
# Check if the input object is in the right format that can be used in the 
# function 

if(class(fit.Model)[1]!="nlme" && class(fit.Model)!="lme" && class(fit.Model)!="lmer")
{
stop("Please specify a object fit by lme or nlme in nlme package or lmer in lme4 package.")
}


###############################################################################
if(!is.null(subset.ids))
{
if(!is.null(pct.rand)) stop("Since \'subset.ids\' was specified, do not also specify \'pct.rand\'.")
if(!is.null(number.rand)) stop("Since \'subset.ids\' was specified, do not also specify \'number.rand\'.")
}

if(!is.null(number.rand))
{
if(!is.null(pct.rand)) stop("Since \'number.rand\' was specified, do not also specify \'pct.rand\'.")
}
###############################################################################

###############################################################################
if(file=="") file <- "vit.fitted"
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



###############################################################################
# This part is for fitted object produced by nlme package

if(class(fit.Model)[1]!="lmer")
{

###############################################################################
if(!requireNamespace("nlme", quietly = TRUE)) stop("The package 'nlme' is needed; please install the package and try again.")
  
all.data <- nlme::augPred(fit.Model)

if(ylab=="") ylab=colnames(all.data)[3]
if(xlab=="") xlab=colnames(all.data)[1]

if(same.scales==TRUE) ylim <- c(min(all.data[3]),max(all.data[3]))
if(same.scales==FALSE) ylim <- NULL

Unique.ID.Values <- unique(as.matrix(all.data[[".group"]]))

if(is.null(number.rand)) 
{
if(is.null(pct.rand)) 
{
pct.rand <- 1
}
if(pct.rand <=1 & pct.rand>0) number.rand <- ceiling(length(Unique.ID.Values)*pct.rand)
if(pct.rand >1 & pct.rand <= 100) number.rand <- ceiling(length(Unique.ID.Values)*pct.rand/100)
}

#if(!is.na(as.numeric(Unique.ID.Values[1])))
#{
#Unique.ID.Values <- as.numeric(Unique.ID.Values)
#}
#
Unique.ID.Values <- sample(Unique.ID.Values, number.rand, replace = FALSE)

if(!is.null(subset.ids))
{

no.ids <- 0
for(i in 1:min(length(Unique.ID.Values),length(subset.ids)))
{
if(sum(Unique.ID.Values==subset.ids[i])==0) no.ids <- cbind(no.ids,subset.ids[i])
}

if(length(no.ids)!=1) stop(paste("id =",paste(no.ids[-1],collapse = ",")," is(are) not in the fitted object."))

Unique.ID.Values <- subset.ids
}

###################################################################################
# N is the number of cases we will plot and in the Quality.of.Fit object

N <- length(Unique.ID.Values)
Residuals <- as.data.frame(cbind(ID=as.character(as.matrix(fit.Model$groups[,1])), 
                                 Y=all.data[all.data[[".type"]]=="original",][,3], 
                                 Fitted=fitted(fit.Model), 
                                 Residual=resid(fit.Model)))

par(mfrow=layout, mar = c(4, 3, 3, 1), mgp = c(2, .5, 0))
plots.per.page <- layout[1]*layout[2]

Quality.of.Fit <- cbind(ID=NA, r2=NA, RMSE=NA)
for(i in 1:N)
{
ids <- Residuals["ID"]==Unique.ID.Values[i]

Y <- as.numeric(as.character(Residuals[ids,"Y"]))
fit <- as.numeric(as.character(Residuals[ids,"Fitted"]))
res <- as.numeric(as.character(Residuals[ids,"Residual"]))

These.data <- all.data[all.data[[".group"]]==Unique.ID.Values[i],]
predicted.data <- These.data[These.data[[".type"]]=="predicted",]
original.data <- These.data[These.data[[".type"]]=="original",]

Quality.of.Fit <- rbind(Quality.of.Fit, c(Unique.ID.Values[i], cor(Y,fit)^2, sqrt(var(res))))
if(i==1) 
{
Quality.of.Fit <- Quality.of.Fit[-1,]
Quality.of.Fit <- t(as.matrix(Quality.of.Fit))
}
ID.i <- Unique.ID.Values[i]
R2.i <- round(as.numeric(Quality.of.Fit[i,2]), 2)
RMSE.i <- round(as.numeric(Quality.of.Fit[i,3]), 2)

plot(original.data[[colnames(These.data)[1]]], original.data[[colnames(These.data)[3]]], 
     main=substitute(italic(R)^2 == R2.i * {" & RMSE" == RMSE.i} * {" (ID" == ID.i} * {")"}, 
     list(ID.i = ID.i, R2.i = R2.i, RMSE.i = RMSE.i)), ylab=ylab, xlab=xlab, ylim=ylim,
     cex.main=.9, cex.lab=.9, ...)


lines(predicted.data[[colnames(These.data)[1]]], predicted.data[[colnames(These.data)[3]]], ...)

}
}



###################################################################################
###################################################################################
###################################################################################
# The following part is for lmer object


if(class(fit.Model)[1]=="lmer")
{

all.data <- attributes(fit.Model)$frame

Unique.ID.Values <- as.matrix(unique(all.data[,dim(all.data)[2]])) 

#dim(all.data)[2]# indicates the location of ID in the data

###############################################################################
if(is.null(number.rand)) 
{
if(is.null(pct.rand)) 
{
pct.rand <- 1
}
if(pct.rand <=1 & pct.rand>0) number.rand <- ceiling(length(Unique.ID.Values)*pct.rand)
if(pct.rand >1 & pct.rand <= 100) number.rand <- ceiling(length(Unique.ID.Values)*pct.rand/100)
}

#if(!is.na(as.numeric(Unique.ID.Values[1])))
#{
#Unique.ID.Values <- as.numeric(Unique.ID.Values)
#}
#
Unique.ID.Values <- sample(Unique.ID.Values, number.rand, replace = FALSE)

if(!is.null(subset.ids))
{

no.ids <- 0
for(i in 1:min(length(Unique.ID.Values),length(subset.ids)))
{
if(sum(Unique.ID.Values==subset.ids[i])==0) no.ids <- cbind(no.ids,subset.ids[i])
}

if(length(no.ids)!=1) stop(paste("id =",paste(no.ids[-1],collapse = ",")," is(are) not in the fitted object."))

Unique.ID.Values <- subset.ids
}

N <- length(Unique.ID.Values)

###############################################################################
coeff <- as.data.frame(coef(fit.Model)[1])

min.max <- c(min(all.data[,2]),max(all.data[,2]))

line.x <- seq(min.max[1],min.max[2],length=100)

par(mfrow=layout, mar = c(4, 3, 3, 1), mgp = c(2, .5, 0))

if(same.scales==TRUE) ylim <- c(min(all.data[,1]),max(all.data[,1]))
if(same.scales==FALSE) ylim <- NULL

if(ylab=="") ylab=names(all.data)[1]
if(xlab=="") xlab=names(all.data)[2]

###############################################################################
Quality.of.Fit <- cbind(ID=NA, r2=NA, RMSE=NA)
y.hat <- NULL
for (i in 1:N)
{
intercept <- coeff[as.matrix(row.names(coeff))==Unique.ID.Values[i],1]
slope <- as.numeric(coeff[as.matrix(row.names(coeff))==Unique.ID.Values[i],-1])

data.x <- all.data[as.matrix(all.data[,dim(all.data)[2]])==Unique.ID.Values[i],2]
data.y <- all.data[as.matrix(all.data[,dim(all.data)[2]])==Unique.ID.Values[i],1]

line.y <- intercept
model.imply.y <- intercept 
for (j in 1)#:length(slope))
{
model.imply.y <- model.imply.y+data.x*slope[j]
line.y <- line.y+line.x*slope[j]
}

residual <- data.y - model.imply.y
R2 <- cor(data.y,model.imply.y )^2
RMSE <- sqrt(var(residual))
ID <- Unique.ID.Values[i]
y.hat <- c(y.hat,model.imply.y)

plot(data.x,data.y, 
     main=substitute(italic(R)^2 == R2 * {" & RMSE" == RMSE} * {" (ID" == ID} * {")"}, 
     list(ID = ID, R2 = round(R2,2), RMSE = round(RMSE,2))), 
     xlim=min.max, ylim=ylim, xlab=xlab,ylab=ylab,
     cex.main=.9, cex.lab=.9, ...)
     
lines(line.x,line.y)

Quality.of.Fit <- rbind(Quality.of.Fit, c(ID, R2, RMSE))
if(i==1) Quality.of.Fit <- Quality.of.Fit[-1,]
}
}

par(mfrow=c(1,1))

###################################################################################
if(save.pdf==TRUE)
{
dev.off()
if(file=="vit.fitted") print(paste("\'vit.fitted.pdf\' file saved at the directory",getwd()))
}

if(save.eps==TRUE)
{
dev.off()
if(file=="vit.fitted") print(paste("\'vit.fitted.eps\' file saved at the directory",getwd()))
}

if(save.jpg==TRUE)
{
dev.off()
if(file=="vit.fitted") print(paste("\'vit.fitted.jpg\' file(s) saved at the directory",getwd()))
}
###################################################################################

Quality.of.Fit <<- cbind(as.data.frame(Quality.of.Fit[,1]),
                         as.data.frame(Quality.of.Fit[,2]),
                         as.data.frame(Quality.of.Fit[,3]))

colnames(Quality.of.Fit) <<- list("ID","R-SQUARE","RMSE")


if(class(fit.Model)[1]!="lmer")
{
observed.data <- all.data[all.data[[".type"]]=="original",][,-4]
fitted.data <- cbind(as.data.frame(observed.data[,2]),as.data.frame(observed.data[,1]),as.data.frame(observed.data[,3]),as.data.frame(fit.Model$fit)[,2])
colnames(fitted.data) <- list(names(fit.Model$groups),xlab,ylab,"fitted")
}

if(class(fit.Model)[1]=="lmer")
{
ID <- as.data.frame(all.data[,3])
fitted.data <- cbind(as.data.frame(all.data[,3]),as.data.frame(all.data[,2]),as.data.frame(all.data[,1]),y.hat)
colnames(fitted.data) <- list(names(all.data)[3],xlab,ylab,"fitted")
}

print("Type \'Quality.of.Fit\' to see the quality of fit.")
print("Type \'fitted.data\' to see the data and fitted data.")

}
