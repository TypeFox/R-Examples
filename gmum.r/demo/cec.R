library(gmum.r)

data(cec.mouse1.spherical)
dataset = cec.mouse1.spherical

# That is the dataset we want to cluster:
plot(dataset, main="Mouse-like dataset")

# Run CEC with default parameters. Set the number of clusters and the dataset.
c <- CEC(k=3, x=dataset)
plot(c)

# Since initial clusterization is random. It may be a good idea to run CEC multiple times and choose the best result.
c <- CEC(k=3, x=dataset, control.nstart=10)
plot(c)

# Better than before, however, we know that clusters are spherical; let's inform CEC about that.
c <- CEC(k=3, x=dataset, control.nstart=10, method.type='spherical')
plot(c)

# Learn details of clustering:
c$centers
c$covMatrix

# Predict cluster which a point would belong to:
predict(c, c(1,1))

# Visualise size and shape of clusters:
plot(c, ellipses=TRUE)

# Try the same with random assignment.
c <- CEC(k=3, x=dataset, control.nstart=10, method.type='spherical', method.init='random')
plot(c)
