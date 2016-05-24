library(ConnMatTools)
data(chile.loco)

num <- prod(dim(chile.loco)) / sum(chile.loco)
betas <- betasVectorDefault(n=num,steps=4)
chile.loco.split <- optimalSplitConnMat(chile.loco,normalize.cols=FALSE,
                                        betas=betas)

# Extra 3rd division
print(paste("Examining split with",names(chile.loco.split$best.splits)[1],
            "subpopulations."))
pops <- subpopsVectorToList(chile.loco.split$subpops[,chile.loco.split$best.splits[[1]]$index])

reduce.loco <- reducedConnMat(pops,chile.loco)

sr <- selfRecruitment(reduce.loco)
lr <- localRetention(reduce.loco)
rlr <- relativeLocalRetention(reduce.loco)
