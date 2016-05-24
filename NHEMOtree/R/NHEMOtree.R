NHEMOtree <-
function(method=c("NHEMO", "NHEMO_Cutoff", "Wrapper"),
                     formula=NULL, data=NULL, 
                     x=NULL, grouping=NULL,
                     CostMatrix, ...){

  # Checking
    method<- match.arg(method)

  # Argument usage
    if ((length(formula)==0 || length(data)==0) && (length(x)==0 || length(grouping)==0)) stop("Incorrect data input!")
    if ((length(formula)!=0 || length(data)!=0) && (length(x)!=0 || length(grouping)!=0)) stop("Incorrect data input!")

  # Transformation of input data
    if (length(formula)>0 && length(data)>0){ # Input by formula and data
        data<- formula_to_dataset(formula, data)
    } else{ # Input by x and grouping
      data<- cbind(grouping, x)
    }

  # Handling of missing data
    dim_data1<- nrow(data)
    data     <- Missing_data(data)  # Deletion of missing data!
    if (dim_data1>nrow(data)) print(paste(dim_data1-nrow(data), "observations deleted due to missing data!"))

  # Transformation of CostMatrix
    if (length(intersect(names(data), CostMatrix[,1]))!=(ncol(data)-1)){ 
        stop("Input data and CostMatrix do not match!")
    }else{
        costs_num<- c()
        for (i in 1:(ncol(data)-1)) costs_num[i]<- which(CostMatrix[,1]==names(data)[i+1])
        CostMatrix<- CostMatrix[costs_num,]
    }

  # Missing data in CostMatrix
    if (ncol(CostMatrix) != 2)                  stop("Incorrect dimension of CostMatrix!")
    if (length(which(is.na(CostMatrix[,2])))>0) stop("Missing data in CostMatrix!")

  # Output
    switch(method,   
      NHEMO        = { # Standard NHEMOtree 
                         result<- NHEMO(data, CostMatrix, ...)
      },

      NHEMO_Cutoff = { # NHEMOtree with local cutoff optimization
                         result<- NHEMO_Cutoff(data, CostMatrix, ...)
      },

      Wrapper      = { # Wrapper approach with NSGA-II and an enclosed classification tree algorithm
                         result<- Wrapper(data, CostMatrix, ...)
    })
    return(result)
}
