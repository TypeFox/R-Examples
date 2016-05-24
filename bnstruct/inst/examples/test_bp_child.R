library("bnstruct")

# child_NA_5000 <- BNDataset(name = "Child")
# child_NA_5000 <- read.dataset(child_NA_5000, header.file="../bnstruct/inst/extdata/Child_data_na_5000.header",
#                               data.file="../bnstruct/inst/extdata/Child_data_na_5000.data",
#                               imputation = TRUE, bootstrap=FALSE, num.boots=3)
# save(child_NA_5000, file="data/child_NA_5000.rda")
# break
# print(mydata)

mydata <- child()

# mydata <- BNDataset("/home/alberto/didattica/tesi/bnstruct/inst/extdata/Child_data_na_5000.data",
#                     "/home/alberto/didattica/tesi/bnstruct/inst/extdata/Child_data_na_5000.header",
#                     starts.from=0)
net <- BN(mydata)
#slot(net, "name") <- "Child"
#slot(net, "num.nodes") <- 20
#slot(net, "node.sizes") <- mydata@node.sizes
#slot(net, "cpts") <- list(NULL)
#slot(net, "variables") <- mydata@variables

# net <- learn.structure(net, mydata, algo="sm", layering= c(1,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5), max.fanin.layers=as.matrix(read.table(header=F,text="
#    0  1  1  1  1
#    0  1  1  1  1
#    0  0  8  7  7
#    0  0  0 14  6
#    0  0  0  0 19")), max.fanin=3)
print(net)
#net <- learn.structure(net, mydata, algo="mmhc", scoring.func = "BIC")
# net <- learn.structure(net, mydata, algo="sm", scoring.func = "BIC",
# layering= c(1,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5), max.fanin.layers=as.matrix(read.table(header=F,text="
#    0  1  1  1  1
#    0  1  1  1  1
#    0  0  8  7  7
#    0  0  0 14  6
#    0  0  0  0 19")), max.fanin=3, bootstrap = FALSE)
#net <- learn.params(net, mydata)

# em(inf.eng, mydata)
#net <- learn.structure(net, mydata, algo="mmhc", scoring.func="BIC")
#net <- learn.params(net, mydata)
print(dag(net))
out <- learn.network(net, mydata, algo = "sem", scoring.func = "BDeu", struct.threshold = 0)
# out <- sem(net, mydata, struct.threshold = 0, algo="mmhc", scoring.func = "BIC")#,
#            layering= c(1,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5), max.fanin.layers=as.matrix(read.table(header=F,text="
#    0  1  1  1  1
#    0  1  1  1  1
#    0  0  8  7  7
#    0  0  0 14  6
#    0  0  0  0 19")), max.fanin=3, bootstrap = FALSE)
# print(out)
# save.to.eps(updated.bn(out$InferenceEngine), filename = "bic.ps")