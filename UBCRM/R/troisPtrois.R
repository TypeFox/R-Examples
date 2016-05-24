troisPtrois <-
function(data = data,lastdose){
idx <- which(data$dose == lastdose)
ndlt <- data$ndlt[idx]
npt <- data$npt[idx]
mtd <- NA
# First study of the dose
if(npt == 3){
nextdose <- ifelse(ndlt == 0, lastdose + 1, ifelse(ndlt > 1, NA, lastdose))
}
# Complementary study of the dose
if(npt == 6){
nextdose <- ifelse(ndlt == 1, lastdose + 1, NA)
}
if(!(nextdose %in% data$dose) | ndlt > 1)
{
nextdose <- NA
mtd <- ifelse(lastdose == 1, NA, ifelse(ndlt > 1, lastdose-1, lastdose))
}
list(nextdose=nextdose,mtd=mtd)
}
