Histogram <- function (vectors, n_classes) 
{
    n_classes = round(360/n_classes)
    abs_his <- rep(0, times = n_classes)
    rel_his <- rep(0, times = n_classes)
    n_vectors = length(vectors)
    for (i in 1:n_vectors) {
        portion = (vectors[i] * n_classes) / 360
        abs_his[portion + 1] = abs_his[portion + 1] + 1
    }
    rel_his = abs_his / n_vectors
    histogram = c(abs_his, rel_his)
    dim(histogram) = c(n_classes, 2)
    return(histogram)
}
