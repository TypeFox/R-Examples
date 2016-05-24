two.stage.model.outputs<-function (B1,catchability,obs,g)
{
#Creation of the dataframe which stores the predicted abundance indices
ai<-data.frame(c(obs$year), c(rep(0,length(obs$year))),c(rep(0,length(obs$year))),c(rep(0,length(obs$year))),c(rep(0,length(obs$year))))
colnames(ai)<-c("year","bts","cgfs","uk","fr")

#loop to compute the predicted abundance indices
for (i in 1:length(obs[,1]))
{
#Prediction for BTS index
ai$bts[i]<-exp(catchability[1])*B1[i]

#Prediction for CGFS index
ai$cgfs[i]<-exp(catchability[2])*B1[i]*exp(-g/4)

#Prediction for UK LPUE
ai$uk[i]<-0.5*exp(catchability[3])*(B1[i]*exp(-g/4)+(B1[i]*exp(-g/2)-obs$L.Q341[i])*exp(-g/4))

#Prediction for FR LPUE
ai$fr[i]<-0.5*exp(catchability[4])*(B1[i] + (B1[i]*exp(-g/2) - obs$L.Q3412[i])*exp(-g/2))
}

#MLH
sigma.bts<-sqrt(mean(log(obs$bts/ai$bts)^2))
sigma.cgfs<-sqrt(mean(log(obs$cgfs/ai$cgfs)^2))
sigma.lpue.uk<-sqrt(mean(log(obs$lpue.uk/ai$uk)^2))
sigma.lpue.fr<-sqrt(mean(log(obs$lpue.fr/ai$fr)^2))

residus.mlh<-data.frame(
c(obs$year[1]:obs$year[length(obs$year)]),
(length(obs$bts)*log(sigma.bts)+1/(2*sigma.bts^2)*sum(log(obs$bts/ai$bts)^2)),
(length(obs$cgfs)*log(sigma.cgfs)+1/(2*sigma.cgfs^2)*sum(log(obs$cgfs/ai$cgfs)^2)),
(length(obs$lpue.uk)*log(sigma.lpue.uk)+1/(2*sigma.lpue.uk^2)*sum(log(obs$lpue.uk/ai$uk)^2)),
(length(obs$lpue.fr)*log(sigma.lpue.fr)+1/(2*sigma.lpue.fr^2)*sum(log(obs$lpue.fr/ai$fr)^2))
)
colnames(residus.mlh)<-c("year","res.bts","res.cgfs","lpue.uk","lpue.fr")

#SSR
residus.ssr<-data.frame(
c(obs$year[1]:obs$year[length(obs$year)]),
(log(obs$bts/ai$bts))^2,
(log(obs$cgfs/ai$cgfs))^2,
(log(obs$lpue.uk/ai$uk))^2,
(log(obs$lpue.fr/ai$fr))^2
)
colnames(residus.ssr)<-c("year","res.bts","res.cgfs","lpue.uk","lpue.fr")


#Residuals (observed - predicted)
residus.raw<-data.frame(
c(obs$year[1]:obs$year[length(obs$year)]),
(obs$bts-ai$bts),
(obs$cgfs-ai$cgfs),
(obs$lpue.uk-ai$uk),
(obs$lpue.fr-ai$fr)
)
colnames(residus.raw)<-c("year","bts","cgfs","uk","fr")

#Standardised residuals
residuals.st<-data.frame(
c(obs$year[1]:obs$year[length(obs$year)]),
(residus.raw$bts-mean(residus.raw$bts))/sd(residus.raw$bts),
(residus.raw$cgfs-mean(residus.raw$cgfs))/sd(residus.raw$cgfs),
(residus.raw$uk-mean(residus.raw$uk))/sd(residus.raw$uk),
(residus.raw$fr-mean(residus.raw$fr))/sd(residus.raw$fr))
colnames(residuals.st)<-c("year","bts","cgfs","lpue.uk","lpue.fr")

#Sum of residuals (to be minimized)
sum.residus<-sum(residus.ssr[,2:5])

#biomass in July (B1, recruitment), in January and in July (B2)
biomass<-data.frame(
c(obs$year[1]:obs$year[length(obs$year)]),
B1,
B1*exp(-g/2),
obs$L.Q3412,
(obs$L.Q3412)/(B1*exp(-g/2)),
(B1*exp(-g/2)-obs$L.Q3412)*exp(-g/2)
)
colnames(biomass)<-c("year","B1","B.jan","landings","exp.rate","B2")

resultat<-list(
obs,
sum.residus,
residus.ssr,
residus.mlh,
residus.raw,
residuals.st,
ai,
biomass
)
names(resultat)<-c("observed", "sum.residuals", "residuals.ssr", "residuals.mlh", "residuals.raw", "residuals.st",  "predicted.ai", "biomass")
return(resultat)
}

######################################################################################

two.stage.model.fit<-function(to.fit, obs.fit, g.fit)
{
resultat.fit<-two.stage.model.outputs(B1=to.fit[1:length(obs.fit$year)], catchability=to.fit[(length(obs.fit$year)+1):(length(obs.fit$year)+4)] ,obs=obs.fit, g=g.fit)
return(resultat.fit$sum.residuals)
}

######################################################################################

delta.glm<-function(input.data)
{
input.data$year<-as.factor(input.data$year)
input.data$fishing.season<-as.factor(input.data$fishing.season)
input.data$rectangle<-as.factor(input.data$rectangle)
input.data$power.class<-as.factor(input.data$power.class)

input.data$factor.month<-"00"
input.data[input.data$month==1,]$factor.month<-"01"
input.data[input.data$month==2,]$factor.month<-"02"
input.data[input.data$month==3,]$factor.month<-"03"
input.data[input.data$month==4,]$factor.month<-"04"
input.data[input.data$month==5,]$factor.month<-"05"
input.data[input.data$month==6,]$factor.month<-"06"
input.data[input.data$month==7,]$factor.month<-"07"
input.data[input.data$month==8,]$factor.month<-"08"
input.data[input.data$month==9,]$factor.month<-"09"
input.data[input.data$month==10,]$factor.month<-"10"
input.data[input.data$month==11,]$factor.month<-"11"
input.data[input.data$month==12,]$factor.month<-"12"

input.data<-data.frame(input.data$year, input.data$fishing.season, input.data$factor.month, input.data$rectangle, input.data$power.class, input.data$lpue)
colnames(input.data)<-c("year", "fishing.season", "month", "rectangle", "power.class", "lpue")

#Creation of presence/absence column called presence
input.data$presence<-1
input.data[input.data$lpue==0,]$presence<-0


#binomial error GLM fitting
binomial.glm<-glm(presence ~ fishing.season + month + rectangle + power.class, family="binomial", data=input.data)

#Summary of the fitting results
binomial.summary<-summary(binomial.glm)

#Residuals
binomial.residuals<-residuals(binomial.glm)

#Fitted values
binomial.fit<-fitted(binomial.glm)

#Gaussian error GLM fitting
positive.input.data<-input.data[input.data$lpue>0,]
gaussian.glm<-glm(log(lpue) ~ fishing.season + month + rectangle + power.class, family="gaussian", data=positive.input.data)

#Summary of the fitting results
gaussian.summary<-summary(gaussian.glm)

#Residuals
gaussian.residuals<-residuals(gaussian.glm)

#Fitted values
gaussian.fit<-fitted(gaussian.glm)

positive.input.data<-input.data[input.data$lpue>0,]
positive.input.data$year<-as.factor(as.character(positive.input.data$year))
positive.input.data$fishing.season<-as.factor(as.character(positive.input.data$fishing.season))
positive.input.data$month<-as.factor(as.character(positive.input.data$month))
positive.input.data$rectangle<-as.factor(as.character(positive.input.data$rectangle))
positive.input.data$power.class<-as.factor(as.character(positive.input.data$power.class))

l.fishing.season<-length(levels(as.factor(positive.input.data$fishing.season)))
l.month<-length(levels(as.factor(positive.input.data$month)))
l.rectangle<-length(levels(as.factor(positive.input.data$rectangle)))
l.power.class<-length(levels(as.factor(positive.input.data$power.class)))
l.total<-l.fishing.season*l.month*l.rectangle*l.power.class

predicted.lpue <- matrix(NA,nrow=l.total,ncol=4)

predicted.lpue[,4] <- rep(levels(positive.input.data$power.class),(l.total/l.power.class))

for(k in 1:l.fishing.season){
for(j in 1:l.rectangle){
for(m in 1:l.month){
start.fishing.season <- k*l.rectangle*l.month*l.power.class - l.rectangle*l.month*l.power.class
end.fishing.season <- k*l.rectangle*l.month*l.power.class
start.rectangle <- start.fishing.season + j*l.month*l.power.class - l.month*l.power.class
end.rectangle <- start.fishing.season + j*l.month*l.power.class
start.month <- start.rectangle + m*l.power.class - l.power.class+1
end.month <- start.rectangle + m*l.power.class
predicted.lpue[start.month:end.month,3] <- rep(levels(positive.input.data$month)[m], l.power.class)
}
predicted.lpue[(start.rectangle+1):end.rectangle,2] <- rep(levels(positive.input.data$rectangle)[j], l.month*l.power.class)
}
predicted.lpue[(start.fishing.season+1):end.fishing.season,1] <- rep(levels(positive.input.data$fishing.season)[k], l.rectangle * l.month* l.power.class)
}

predicted.lpue<-as.data.frame(predicted.lpue)
colnames(predicted.lpue)<-c("fishing.season", "rectangle", "month", "power.class")

binomial.glm.prediction <- predict.glm(binomial.glm, predicted.lpue[,1:4], type="r", se.fit=T)
gaussian.glm.prediction <- predict.glm(gaussian.glm, predicted.lpue[,1:4], type="r", se.fit=T)

binomial.glm.prediction.fit<-binomial.glm.prediction$fit
gaussian.glm.prediction.fit<-gaussian.glm.prediction$fit

gaussian.glm.prediction.var<-(gaussian.glm.prediction$se.fit)^2

gaussian.prediction<-exp(gaussian.glm.prediction.fit+gaussian.glm.prediction.var/2)

predicted.lpue$st.lpue<-binomial.glm.prediction.fit*gaussian.prediction

delta.outputs<-list(
binomial.glm,
binomial.summary,
binomial.residuals,
binomial.fit,
gaussian.glm,
gaussian.summary,
gaussian.residuals,
gaussian.fit,
predicted.lpue
)

names(delta.outputs)<-c("binomial.glm",
"binomial.summary",
"binomial.residuals",
"binomial.fit",
"gaussian.glm",
"gaussian.summary",
"gaussian.residuals",
"gaussian.fit",
"predicted.lpue"
)

return(delta.outputs)
}

