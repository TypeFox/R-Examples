library("bnstruct")

# asia_10000 <- BNDataset(name = "Asia")
# asia_10000 <- read.dataset(asia_10000, "../bnstruct/inst/extdata/asia_10000.header", "../bnstruct/inst/extdata/asia_10000.data", imputation = FALSE, bootstrap=FALSE)
# save(asia_10000, file="data/asia_10000.rda")
# break
mydata <- asia()
#mydata@name <- "Asia"
# print(mydata)

#other <- BNDataset(data = mydata@raw.data, variables = variables(mydata), node.sizes = mydata@node.sizes, discreteness = rep('d',8))
#print(other)
# mydata <- BNDataset("/home/alberto/didattica/tesi/bnstruct/inst/extdata/asia_10000.data",
#                     "/home/alberto/didattica/tesi/bnstruct/inst/extdata/asia_10000.header",
#                     starts.from=0)
net <- learn.network(mydata) #, algo = "sm", scoring.func = "BIC")
# set.name(net) <- "Asia"
# set.num.nodes(net) <- 8
# set.node.sizes(net) <- c(2,2,2,2,2,2,2,2)
# slot(net, "cpts") <- list(NULL)
# slot(net, "variables") <- c("Asia", "Tubercolosys", "Smoke", "LungCancer", "Bronchitis", "Either", "X-ray", "Dyspnea")
# set.dag(net) <- matrix(c(0, 1, 0, 0, 0, 0, 0, 0,
#                          0, 0, 0, 0, 0, 1, 0, 0,
#                          0, 0, 0, 1, 1, 0, 0, 0,
#                          0, 0, 0, 0, 0, 1, 0, 0,
#                          0, 0, 0, 0, 0, 0, 0, 1,
#                          0, 0, 0, 0, 0, 0, 1, 1,
#                          0, 0, 0, 0, 0, 0, 0, 0,
#                          0, 0, 0, 0, 0, 0, 0, 0), nrow = 8, ncol = 8, byrow=TRUE)
# 
# 
# dim.vars <- NULL
# dim.vars[[1]] <- as.list(c(1))
# dim.vars[[2]] <- as.list(c(2,1))
# dim.vars[[3]] <- as.list(c(3))
# dim.vars[[4]] <- as.list(c(4,3))
# dim.vars[[5]] <- as.list(c(5,3))
# dim.vars[[6]] <- as.list(c(6,2,4))
# dim.vars[[7]] <- as.list(c(7,6))
# dim.vars[[8]] <- as.list(c(8,5,6))
# 
# cpts <- NULL
# cpts[[1]] <- array(c(0.01, 0.99), dim=c(2), dimnames=list("Asia" = c(1,2)))
# cpts[[2]] <- array(c(0.05, 0.95, 0.01, 0.99), dim=c(2,2), dimnames=list("Tubercolosys" = c(1,2), "Asia" = c(1,2)))
# cpts[[3]] <- array(c(0.5, 0.5), dim=c(2), dimnames=list("Smoke" = c(1,2)))
# cpts[[4]] <- array(c(0.1, 0.9, 0.01, 0.99), dim=c(2,2), dimnames=list("LungCancer" = c(1,2), "Smoke" = c(1,2)))
# cpts[[5]] <- array(c(0.6, 0.4, 0.3, 0.7), dim=c(2,2), dimnames=list("Bronchitis" = c(1,2), "Smoke" = c(1,2)))
# cpts[[6]] <- array(c(1,0,1,0,1,0,0,1), dim=c(2,2,2), dimnames=list("Either" = c(1,2), "Tubercolosys" = c(1,2), "LungCancer" = c(1,2)))
# cpts[[7]] <- array(c(0.98, 0.02, 0.05, 0.95), dim=c(2,2), dimnames=list("X-ray" = c(1,2), "Either" = c(1,2)))
# cpts[[8]] <- array(c(0.9, 0.1, 0.8, 0.2, 0.7, 0.3, 0.1, 0.9), dim=c(2,2,2), dimnames=list("Dyspnea" = c(1,2), "Bronchitis" = c(1,2), "Either" = c(1,2)))
# # 
# slot(net, "cpts") <- cpts

#net <- learn.structure(net, mydata, algo="mmhc")
#net <- learn.params(net, mydata)
 
# print(net)
# unlist(lapply(net@cpts, sum))
# readLines(file("stdin"),1)


print(net)

inf.eng <- InferenceEngine(net)
#inf.eng <- build.junction.tree(inf.eng, net@dag)

print(inf.eng)

observations(inf.eng) <- list(c("Asia", "X-ray", "Dyspnea"), c(1,1,1))
inf.eng <- belief.propagation(inf.eng)

print("-------------------------------------------------------------------------------------------")

marginals(inf.eng)
get.most.probable.values(inf.eng)
