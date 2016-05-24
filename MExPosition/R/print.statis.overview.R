print.statis.overview <- function (x,...) {

  res.statis.overview <- x
  if (!inherits(res.statis.overview, "statis.overview")) stop ("no convenient data")
  cat("**Overview of Results for STATIS **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(8, 2), list(1:8, c("Name", "Description")))
  
  res[1,] <- c("$data","Raw data matrix")
  res[2,] <- c("$groupmatrix","Group matrix representing the different tables")
  res[3,] <- c("$preprocess.data","Preprocessed data matrix")
  res[4,] <- c("$num.groups","Number of groups")
  res[5,] <- c("$num.obs","Number of observations")
  res[6,] <- c("$row.preprocess","Row preprocessing option selected")
  res[7,] <- c("$column.preprocess","Column preprocessing option selected")
  res[8,] <- c("$table.preprocess","Table preprocessing option selected")
  
  print(res)
}
