require("GMD") # load library
data(cage)     # load data

## measure pairwise distance
x <- gmdp(cage[["Pfkfb3 (T02R00AEC2D8)"]],cage[["Csf1 (T03R0672174D)"]])
print(x)                     # print a brief version by default
print(x, mode="full")  # print a full version by default

## show alignment
plot(x,labels=c("Pfkfb3","Csf1"),beside=FALSE)

## show another alignment
plot(gmdp(cage[["Hig1 (T09R0743763C)"]],cage[["Cd72 (T04R028B8BC9)"]]),
     labels=c("Hig1 (T09R0743763C)","Cd72 (T04R028B8BC9)"),
     beside=FALSE)

## construct a distance matrix and visualize it
short.labels <- gsub("(.+) \\(.+","\\1",names(cage)) # get short labels
x <- gmdm(cage[1:6],labels=short.labels[1:6])
plot(x)
