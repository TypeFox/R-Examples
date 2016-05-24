################################
## Example for STdata/STmodel ##
################################
##load data
data(mesa.model)

##standard plot
qqnorm(mesa.model)
##add a line, and group (and colour) by AQS/FIXED
par(mfrow=c(2,2))
obs.type <- mesa.model$locations$type[match(mesa.model$obs$ID,
                                            mesa.model$locations$ID)]
qqnorm(mesa.model, line=1, group=obs.type, col=obs.type)

##colour code by season and split by type
##First create a vector dividing data into four seasons
I.season <- as.factor(as.POSIXlt(mesa.model$obs$date)$mon+1)
levels(I.season) <- c(rep("Winter",2), rep("Spring",3), 
                      rep("Summer",3), rep("Fall",3), "Winter") 

par(mfrow=c(2,2))
qqnorm(mesa.model, line=1, col=I.season, group=obs.type)
legend("bottomright", legend=as.character(levels(I.season)),
       pch=1, col=1:nlevels(I.season))

###############################
## Example for predCVSTmodel ##
###############################
##load data
data(pred.cv.mesa)

##standard plot
par(mfrow=c(1,1))
qqnorm(pred.cv.mesa, line=3)
##or for the normalised residuals
qqnorm(pred.cv.mesa, line=3, norm=TRUE)

##add a line, and group by AQS/FIXED
par(mfrow=c(2,2))
qqnorm(pred.cv.mesa, line=1, group=obs.type)

##and for normalised residuals, colour-coded by season
par(mfrow=c(2,2))
qqnorm(pred.cv.mesa, line=2, norm=TRUE,
       group=obs.type, col=I.season)
legend("bottomright", legend=as.character(levels(I.season)),
       pch=1, col=1:nlevels(I.season))
