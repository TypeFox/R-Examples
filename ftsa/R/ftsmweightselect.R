ftsmweightselect <-
function(data, ncomp=6, ntrainyear, errorcriterion=c("mae","mse","mape"))
{
errorcriterion = match.arg(errorcriterion)
trainyear = max(data$time)-ntrainyear
testdat = extract(data, direction = "time", timeorder = (max(data$time)-ntrainyear+1):max(data$time))$y
if(errorcriterion == "mae")
{
optifunction = function(beta)
{
traindata = matrix(,length(data$x),ntrainyear)
for(i in 1:ntrainyear)
{
traindata[,i] =  forecast(ftsm(extract(data, direction = "time", 
  timeorder = data$time[1]:trainyear), order = ncomp,
  weight = TRUE, beta = beta), h = 1)$mean$y
trainyear = trainyear + 1
}
return(mae(traindata, testdat))
}
weit = optimize(optifunction,c(0,1))$minimum
}
if(errorcriterion == "mse")
{
optifunction = function(beta)
{
traindata = matrix(,length(data$x),ntrainyear)
for(i in 1:ntrainyear)
{
traindata[,i] =  forecast(ftsm(extract(data, direction = "time", 
  timeorder = data$time[1]:trainyear), order = ncomp,
  weight = TRUE, beta = beta), h = 1)$mean$y
trainyear = trainyear + 1
}
return(mse(traindata, testdat))
}
weit = optimize(optifunction,c(0,1))$minimum
}
if(errorcriterion == "mape")
{
optifunction = function(beta)
{
traindata = matrix(,length(data$x),ntrainyear)
for(i in 1:ntrainyear)
{
traindata[,i] =  forecast(ftsm(extract(data, direction = "time", 
  timeorder = data$time[1]:trainyear), order = ncomp,
  weight = TRUE, beta = beta), h = 1)$mean$y
trainyear = trainyear + 1
}
return(mape(traindata, testdat))
}
weit = optimize(optifunction,c(0,1))$minimum
}
return(weit)
}

