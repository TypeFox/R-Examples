LoadData <- function(FileName, Type = 2)
{
	data_=ReadFromFile(FileName);
	
    if (length(data_) > 1) {
      if (CorrectType(Type, data_) == FALSE){
        switch (Type, print("Error, the file is not of Cartesian coordinates type (type 1)."),
				      print("Error, the file is not of Incremental data type (type 2)."),
            		  print("Error, the file is not of Polar coordinates type (type 3)."),
                );
      } else {  
        
      if (Type==1){
        increm=ToCalculateIncr(data_);
        polar_vectors=VectorsToPolar(increm); 
        rectangular_vectors=increm;
      }
	  
	  if (Type==2){
          polar_vectors=VectorsToPolar(data_);	
          rectangular_vectors=data_;
	  }
      
      if(Type==3){
          rectangular_vectors=VectorsToRectangular(data_);
          polar_vectors=data_;
	  }
        
        num_data = dim(data_);
        res = matrix(nrow = num_data[1], ncol = 9); 
        res[,1] = polar_vectors[,1];
        res[,2] = polar_vectors[,2];
        res[,3] = rectangular_vectors[,1];
        res[,4] = rectangular_vectors[,2];
        
        if(Type==1)  {        
          res[,5]=data_[,1];
          res[,6]=data_[,2];
          res[,7]=data_[,3];
          res[,8]=data_[,4];
        }
        return(res)
      }
	}    
  }