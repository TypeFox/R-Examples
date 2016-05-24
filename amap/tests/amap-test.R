
library(amap)

 set.seed(1234)

data(USArrests)

  METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
               "binary","pearson","correlation","spearman","kendall",
               "abspearson","abscorrelation")
  METHODSLINKS <- c("ward", "single", "complete", "average", "mcquitty", 
                    "median", "centroid","centroid2")



for (mymethod in METHODS) {		    
    d = Dist(USArrests, method = mymethod)
    print(d)
    k  = Kmeans(USArrests, centers = 4, method = mymethod)
    print(k)
    for (mylink in METHODSLINKS)
    {
	cat(mylink)
	cat(mymethod)
	hc <- hcluster(USArrests,link = mylink, method =  mymethod, nbproc=4)
	print(hc)
    }
}

hc <- hcluster(USArrests, nbproc=1)
print(hc)
    





KERNELS = c("gaussien", "quartic", "triweight", "epanechikov" , 
"cosinus", "uniform")

for(myKernel in KERNELS) {
  myacp = acprob(USArrests, kernel = myKernel);
  print(myacp)
} 

x = matrix (1:12, 3, 4)
ncol = as.integer(dim(x)[2])
nrow = as.integer(dim(x)[1])
rien = .C("checkMatrix",as.double(x), nrow, ncol, PACKAGE="amap")

m = matrix(1:16,4,4)
d = as.dist(m)
rien = .C("checkMatrixTriangle",as.double(d),as.integer(4), as.integer(0))

d=c(1,2,3,4,6,7,8,11,12,16)
rien = .C("checkMatrixTriangle",as.double(d),as.integer(4), as.integer(1))



