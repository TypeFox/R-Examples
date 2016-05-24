print.doact.statis.overview <- function(x,...){
	
	res.doact.statis.overview <- x
	if (!inherits(res.doact.statis.overview, "doact.statis.overview")) stop ("no convenient data")
  cat("**Oveview of the Analysis via DO-ACT **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(12, 2), list(1:12, c("Name", "Description")))
	
	res[1,] <- c("$data1", "Dataset 1")
	res[2,] <- c("$data2", "Dataset 2")
	res[3,] <- c("$column.desig.1","Column Design for Dataset 1")
	res[4,] <- c("$column.design.2","Column Design for Dataset 2")
	res[5,] <- c("$row.preprocess.data1","Row preprocessing option for Dataset 1")
	res[6,] <- c("$column.preprocess.data1","Colum preprocessing option for Dataset 1")
	res[7,] <- c("$table.preprocess.data1","Table preprocessing option for Dataset 1")
	res[8,] <- c("$row.preprocess.data2","Row preprocessing option for Dataset 2")
	res[9,] <- c("$column.preprocess.data2","Column preprocessing option for Dataset 2")
	res[10,] <- c("$table.preprocess.data2","Table preprocessing option for Dataset 2")
	res[11,] <- c("$num.groups.1","Number of groups in Dataset 1")
	res[12,] <- c("$num.groups2","Number of groups in Dataset 2")
	
	print(res)
}