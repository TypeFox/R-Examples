test.search <- function () {
  for (model in 1:8) {
    result <- searchCrossOverDesign(s=9, p=5, v=4, model=4, eff.factor=1, n=c(25,1))
  }
  for (model in 1:8) {
    v <- 4
    result <- searchCrossOverDesign(s=9, p=5, v=v, model=4, eff.factor=1, n=c(25,1), balance.p=TRUE)
    for (i in 1:v) {
      v.count <- apply(getDesign(result), 1, function(x) {sum(x==i)})
      checkTrue(max(v.count)-min(v.count)<1.5)
    } 
    result <- searchCrossOverDesign(s=9, p=5, v=v, model=4, eff.factor=1, n=c(25,1), balance.s=TRUE)
    for (i in 1:v) {
      v.count <- apply(getDesign(result), 2, function(x) {sum(x==i)})
      checkTrue(max(v.count)-min(v.count)<1.5)
    } 
  }  
}

test.random.matrix.generation <- function() {
  v <- 4
  for (j in 1:5) {
    design <- randomDesign(s=9, p=5, v=4,  balance.p=TRUE, model=1)
    for (i in 1:v) {
      v.count <- apply(design, 1, function(x) {sum(x==i)})
      checkTrue(max(v.count)-min(v.count)<1.5)
    } 
    design <- randomDesign(s=9, p=5, v=4,  balance.s=TRUE, model=1) 
    for (i in 1:v) {
      v.count <- apply(design, 2, function(x) {sum(x==i)})
      checkTrue(max(v.count)-min(v.count)<1.5)
    } 
  }
}

test.strangeDesignInputs <- function() {
  s <- 4 # number of sequences
  p <- 4 # number of periods
  v <- 4 # number of treatments
  
  D <- rbind(c("A","B","C","D"),
             c("B","C","D","A"),
             c("C","D","A","B"),
             c("D","A","B","C"))
  
  D <- matrix(as.numeric(as.factor(D)), dim(D)[1])  
  
  myInv <- ginv(rcd(D, v, model=1))
  
}