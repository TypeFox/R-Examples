### model fake parameters
model.fake.par <- function(fake.model=c("uninformative","average","slight","extreme")) {

  fake.model <- match.arg(fake.model)
  
  if (fake.model=="uninformative") {
    gam <- c(1,1); del <- c(1,1)
  }
  if (fake.model=="average") {
    gam <- c(3,3); del <- c(3,3)
  }
  if (fake.model=="slight") {
    gam <- c(1.5,4); del <- c(4,1.5)
  }
  if (fake.model=="extreme") {
    gam <- c(4,1.5); del <- c(1.5,4)
  }
  
  return(list(gam=gam,del=del))
}