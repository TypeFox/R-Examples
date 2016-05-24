## ---- echo=TRUE, results='asis', eval=TRUE-------------------------------
library(lasvmR)

# generate two gaussian clusters (very easy example)
qx = rnorm(100, mean = -3, sd = 1) - 1
qy = rnorm(100, mean = -3, sd = 1) - 1
px = rnorm(100, mean = 3, sd = 1) + 1
py = rnorm(100, mean = 3, sd = 1) + 1
data = rbind( cbind(px, py), cbind(qx, qy) )
label = sign (data[,1])
model = lasvmTrain (x = data, y = label, gamma = 1.0, cost = 1.0, epochs = 33, epsilon = 0.01)
result = lasvmPredict (data, model)

# compute error
error = sum(abs(label - result$predictions))/length(label)*100
cat ("Error rate is ", error)

