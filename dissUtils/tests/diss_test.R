require(dissUtils);

mat <- matrix(c(1,0,1,0.5,0,1,0,0,0.5,0.5,1,0,1,0,0.5),ncol=3);

approaches <- c("braycurtis",
                "canberra",
                "chebyshev",
                "covariance",
                "euclidean",
                "equality",
                "hellinger",
                "jaccard",
                "mahalanobis",
                "manhattan",
                "minkowski",
                "pearson",
                "procrustes");

for(a in approaches){
    cat(a,"\n");
    print(diss(mat, method = a), digits = 3);
}
