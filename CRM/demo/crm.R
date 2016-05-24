
#The table 1 in O'Quingley et al.'s paper, page 40
#This example is used to illustrate how the program is used to find
#the MTD and the updated parameter 
 
target <- 0.2
prior <- c(0.05,0.1,0.2,0.3,0.5,0.7)
x <- c(3,4,4,3,3,2,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1)
y <- c(0,0,1,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,1)
ptdata <- cbind(x,y)
for(i in 1:25){
  if(i == 1){
    cat(1,1,3,0,"\n")
  }
  res <- crm(target,prior,ptdata[1:i,],model=1,a0=1)
  if(i < 25){
    cat(i+1,res$a,res$MTD,ptdata[i+1,2],"\n")
  }else {
    cat(i+1,res$a,res$MTD,"\n")
  }
}

#the proposed MTD is
res$MTD
