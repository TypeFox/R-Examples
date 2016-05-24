##load the data
data(mesa.model)

##create a matrix with the CV-schemes
I.cv <- createCV(mesa.model, groups=10)

##number of observations in each CV-group
table(I.cv)

##Which sites belong to which groups?
ID.cv <- sapply(split(mesa.model$obs$ID, I.cv),unique)
print(ID.cv)

##Note that the sites with distance 0.084<min.dist 
##are grouped together (in group 10).
mesa.model$D.beta[ID.cv[[10]], ID.cv[[10]]]

##Find out which location belongs to which cv group
I.col <- apply(sapply(ID.cv,function(x) mesa.model$locations$ID
               %in% x), 1, function(x) if(sum(x)==1) which(x) else 0)
names(I.col) <- mesa.model$locations$ID
print(I.col)

##Plot the locations, colour coded by CV-grouping
plot(mesa.model$locations$long, mesa.model$locations$lat,
     pch=23+floor(I.col/max(I.col)+.5), bg=I.col, 
     xlab="Longitude", ylab="Latitude")

###############################################################
## Using matrix representation of cross-validation structure ##
###############################################################

##create a matrix with the CV-schemes
I.cv <- createCV(mesa.model, groups=10, Icv.vector=FALSE)

##number of observations in each CV-group
colSums(I.cv)

##Which sites belong to which groups?
ID.cv <- apply(I.cv, 2, function(x){ unique(mesa.model$obs$ID[x]) })

##and then as above...
