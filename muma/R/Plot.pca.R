Plot.pca <-
function(pcx,pcy,scaling, test.outlier=TRUE) {
 Plot.pca.score(pcx,pcy,scaling)
 Plot.pca.loading(pcx,pcy,scaling)
if (test.outlier) {outlier(pcx, pcy, scaling)}
}
