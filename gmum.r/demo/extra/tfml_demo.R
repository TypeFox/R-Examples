library(gmum.r)

# --- CEC

# Load a data set
data(cec_mouse_1_spherical)
mouse <- input
par(mfrow=c(1,1))

# See what are we dealing with
plot(mouse, main="Mouse dataset")

# Run CEC on it
c <- CEC(k=3, x=mouse)
plot(c)
title(main="CEC on mouse")

# Since initial clusterization is random. 
# It may be a good idea to run cec multiple times 
# and choose the best result.
c <- CEC(k=3, x=mouse, control.nstart=10)
plot(c)
title(main="Best result from 10 CEC runs")

# Better than before, however, we know that clusters 
# are spherical; let's inform cec about that.
c <- CEC(k=3, x=mouse, method.type='spherical')
plot(c)
title(main="CEC with predestined cluster type")

# You can learn details of clustering like this
centers(c)
cov(c)

# You can visualise size and shape of clusters
plot(c, ellipses=TRUE)
title(main="CEC cluster shapes")

# Unfair k-means comparison
library(stats)
cec <- CEC(k=3, x=mouse, control.nstart=10)
km <- kmeans(centers=3, x=mouse, nstart=10)

par(mfrow=c(1,2))
plot(mouse, col=km$cluster)
title(main="k-means")
plot(cec)
title(main="CEC")

# --- SVM

# Data sets
data(svm_two_circles_dataset)
data(svm_two_ellipsoids_dataset)

# Let's run 2 SVMs on first data set, 
# one normal, the other with 2e preprocessing
svm <- SVM(V3~., 
           svm.twoellipsoids.dataset,
           verbosity = 0);

esvm <- SVM(V3~.,
            svm.twoellipsoids.dataset,
            prep="2e",
            verbosity = 0);

# Plot the results
p1 <- plot(svm) +
      scale_x_continuous(limits=c(-10, 10)) + 
      scale_y_continuous(limits=c(-10, 10)) + 
      ggtitle("SVM on two elipsoids")
p2 <- plot(esvm) +
      scale_x_continuous(limits=c(-10, 10)) +
      scale_y_continuous(limits=c(-10, 10)) +
      ggtitle("2eSVM on two elipsoids")
multiplot(p1,p2)

# Let's try the same thing on the other data set
svm <- SVM(V3~.,
           svm.twocircles.dataset,
           verbosity = 0);
esvm <- SVM(V3~.,
            svm.twocircles.dataset, 
            prep      = "2e",
            verbosity = 0);

# Plot the results
p1 <- plot(svm) +
  scale_x_continuous(limits=c(-5, 10)) + 
  scale_y_continuous(limits=c(-5, 5)) + 
  ggtitle("SVM on two elipsoids")
p2 <- plot(esvm) +
  scale_x_continuous(limits=c(-5, 10)) +
  scale_y_continuous(limits=c(-5, 5)) +
  ggtitle("2eSVM on two elipsoids")
multiplot(p1,p2)

# ---

# Let's say we changed our mind
svm <- SVM(V3~., 
           svm.twocircles.dataset, 
           lib="svmlight",
           verbosity=0);
plot(svm)
svm$call

# ---

# Load some data and weights
load("demo/extra/df.RDa")
load("demo/extra/w.RDa")
formula <- Y~.
const_weights <- rep(1,nrow(df))

# Fit two SVMs, one with weighted examples
svm <- SVM(formula, 
           df, 
           lib       = "svmlight", 
           kernel    = "rbf", 
           C         = 10, 
           gamma     = 30, 
           verbosity = 0);

wsvm <- SVM(formula, 
            df, 
            lib             = "svmlight", 
            kernel          = "rbf", 
            example.weights = weights, 
            C               = 10, 
            gamma           = 30, 
            verbosity       = 0);

# Plot
p1 <-plot(svm)
p2 <- plot(wsvm)

