print.RWBP <-
function(x,...){
cat(" A Random Walk on Bipartite Graph spatial outlier detection model was built: \n",
"---------------------------------------------------------------------------- \n\n",
"neighberhood size = ", x$nn_k, "\n","initial clusters amount = ", x$k, "\n",
"each process increases clusters amount by ", x$clusters.stepSize, " more clusters\n",
"clusters iterations amount = ", x$h, "\n", "alfa = ", x$alfa, "\n",
"dumping factor = ", x$c, "\n", "valid rows = ", x$n," out of ", x$n.orig," input rows (records with empty values were removed)", "\n")
cat("\n","a bipartite graph was built: \n")
print(x$igraph)
cat("\n","outlier scores: ","\n")
print(x$OutScore)
}
