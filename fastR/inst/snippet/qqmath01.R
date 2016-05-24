p <- qqmath(~Sepal.Length|Species,data=iris)
print(p)
set.seed(1)                      # use fixed random seed
p <- qqmath(~rnorm(150)|factor(rep(1:3,each=50)))
print(p)
