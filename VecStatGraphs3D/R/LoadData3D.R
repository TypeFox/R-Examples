LoadData3D <- function (FileName, Type = 2) 
{
    data_ = ReadFromFile3D(FileName)
	
    if (length(data_) > 1) {
        if (CorrectType3D(Type, data_) == FALSE) {
            switch(Type, print("Error, expected file is in the format X origin, Y origin, Z origin, X end, Y end, Z end"), 
                print("Error, expected file is in the format X increment, Y increment, Z increment"), 
                print("Error, expected file is in the format module, colatitude, longitude"))
        }
        else {
            if (Type == 1) {
				axis_incr = ToCalculateIncr3D(data_)
				polar_vectors = VectorsToPolar3D(axis_incr)
				rectangular_vectors = axis_incr
            }
            if (Type == 2) {
                polar_vectors = VectorsToPolar3D(data_)
                rectangular_vectors = data_
            }
            if (Type == 3) {
	            polar_vectors = data_
                rectangular_vectors = VectorsToRectangular3D(polar_vectors)
            }
			
			x <- rectangular_vectors[, 1]
            y <- rectangular_vectors[, 2]
            z <- rectangular_vectors[, 3]
            num_data = dim(data_)
            res = matrix(nrow = num_data[1], ncol = 12)
            res[, 1] = 1
            res[, 2] = polar_vectors[, 1]
            res[, 3] = polar_vectors[, 2]
            res[, 4] = rectangular_vectors[, 1]
            res[, 5] = rectangular_vectors[, 2]
            res[, 6] = rectangular_vectors[, 3]
            if (Type == 1) {
                res[, 7] = data_[, 1]
                res[, 8] = data_[, 2]
                res[, 9] = data_[, 3]
                res[, 10] = data_[, 4]
                res[, 11] = data_[, 5]
                res[, 12] = data_[, 6]
            }
            return(res)
        }
    }
}
