sim3p3 <-
function(truerate, seed=NULL){
if (!is.null(seed)) {set.seed(seed)}
# prob = probability vector
data <- CreData(length(truerate))
nextdose <- 1
while (nextdose %in% data$dose) {
lastdose <- nextdose
ndlt <- sum(rbinom(3, 1, truerate[lastdose]))
data <- updata(data, lastdose, 3, ndlt)
nextdose <- troisPtrois(data, lastdose)$nextdose
}
mtd <- troisPtrois(data, lastdose)$mtd
list(data=data,mtd=mtd,lastdose=lastdose)
}
