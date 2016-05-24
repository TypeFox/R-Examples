# This function allows manual entry of the data and returns a matrix sorted by the X values
# in ascending order
fillData <- function(){
	# Create a data frame and then open it to manual edit
	mydata <- data.frame(X=numeric(0),Y=numeric(0),Weight=numeric(0))
	mydata <- edit(mydata)
	# Sort the data by the X value
	mydata<-mydata[order(mydata[1]),]
	return (mydata)
}

# This funciton receives a matrix and retruns the number of lines in it. The number of
# lines is derived from the number of X values present in the matrix
size<-function(mat){
	i<-1
	while(!is.na(mat[i,1])&!is.null(mat[i,1])){
		i=i+1
	}
	return(i-1)
}

# This function checks if weights were assigned to all the data pairs (x,y)
# if weights were assigned, then the function replaces this weight with the normalized
# weight. If weights weren't provided, equal weights are provided to all pairs
calcWeights<-function(mat){
	n<-size(mat)
	zeroWeight<-0
	sumWeight<-0
	for (i in 1:n){
		# If this observation has no weight specified, update the flag zeroWeight
		if (is.null(mat[i,3])||is.na(mat[i,3])||mat[i,3]==0){
			zeroWeight<-1
		}
		sumWeight<-sumWeight+mat[i,3]
	}
	# If at least one observation pair doesn't have a weight specified, ignore
	# all the weights and assume all the observations are of equal weight
	if (zeroWeight==1){
		mat[3]=1/n
		return (mat)
	}
	# Else calculate the relative weight of each observation
	for (i in 1:n){
		mat[i,3]<-mat[i,3]/sumWeight
	}
	return(mat)
}

# This function calculates the estimate for the cumulative distribution function of each of
# the X values and stores this value in the fourth column of the matrix
calcFX<-function(mat){
	mat[1,4]=mat[1,3]
	n<-size(mat)
	for (i in 2:n){
		# The F(X) of the current observation is equal to the F(X) of the previous
		# observation plus the weight of the current observation
		mat[i,4]=mat[i-1,4]+mat[i,3]
	}
	names(mat)[4]<-"FX"
	return (mat)
}

# This function calculates the cumulative average value for each Y value and stores this 
# values in the 5th column of the matrix
calcFY<-function(mat){
	mat[1,5]=mat[1,2]*mat[1,3]
	n<-size(mat)
	for (i in 2:n){
		# The F(Y) value of the current observation is equal to the F(Y) calculated
		# for the previous observation plus the current observation Y value
		# multiplied by it's weight
		mat[i,5]=mat[i-1,5]+mat[i,2]*mat[i,3]
	}
	names(mat)[5]<-"FY"
	return (mat)
}

# This function calculates the value of the LOI for each of the matrix rows and stores
# it in the 6th column of the matrix
calcLOI<-function(mat){
	n<-size(mat)
	for (i in 1:n){
		# The fitting LOI value for each observation is the average Y value of the
		# entire sample multiplied by the corresponding F(X) value
		mat[i,6]=mat[i,4]*mat[n,5]
	}
	names(mat)[6]="LOI"
	return (mat)
}

# This function calculates the LMA (LOI minus ACC) for each of the matrix lines and stores
# this value in the 7th column of the matrix
calcLMA<-function(mat){
	n<-size(mat)
	for (i in 1:n){
		mat[i,7]=mat[i,6]-mat[i,5]
	}
	names(mat)[7]="LMA"
	return (mat)
}

# This function plots the ACC/LOE and the LMA graphs base on the data stored in the matrix
# passed as the mat argument
plotGraphs<-function(mat){
	# Adding an all zeroes rows to the matrix so that the curves will start from
	# the axis origin
	mat[size(mat)+1,4]=0
	mat[size(mat)+1,5]=0
	mat[size(mat)+1,6]=0
	mat[size(mat)+1,7]=0
	# Then sorting the matrix according to the F(X) values (this is done in order
	# to move the zeroes row added above to the beginning of the matrix)
	mat<-mat[order(mat[4]),]
	# Then the matrix is transformed making it easier to work with in the plotting
	# functions
	trans<-t(mat)
	# Saving the original plotting settings of the system (colors and etc)
	originalPar<-par(no.readonly = TRUE)
	# Defining custom plotting settings for the ACC vs LOI graph
	par(lwd=2)
	par(col="black")
	# Plotting the ACC vs LOI graph
	plot(trans[4,],trans[5,],type="n",main="ACC",xlab="F(x)",ylab="Accumulating Y Mean")
	# Defining the line colors and data source (the curves) of the ACC vs LOI graph
	par(col="blue")
	lines(trans[4,],trans[5,],type="l")
	par(col="red")
	lines(trans[4,],trans[6,],type="l")
	# Defining the legend
	legend("bottomright",c("ACC","LOE"),col=c("blue","red"),lwd=2,bty="n",text.col="black")
	# Restoring the default plotting settings
	par(originalPar)
	# Opening a new window to plot the LMA graph
	dev.new()
	# Defining custom plotting settings for the LMA graph
	par(lwd=2)
	par(col="black")
	par(xaxs="i")
	# Plotting the LMA graph
	plot(trans[4,],trans[7,],type="n",main="LMA",xlab="F(x)",ylab="LOE minus ACC")
	# Defining the line colors and data source (the curves) of the LMA graph
	par(col="black",lwd=1)
	lines(c(0,1),c(0,0),type="l")
	par(col="blue",lwd=2)
	lines(trans[4,],trans[7,],type="l")
	# Restoring the default plotting settings
	par(originalPar)
}

# This function opens a csv file and copies it contents into a matrix. This function assumes
# that the first column in the file contain the X values, the second column the Y values and
# the 3rd column the relevant weight of each line. The head parameter defines by TRUE/FALSE
# if the file contains a header or not
fillCSVData<-function(str,head){
	# Read data from the file to a matrix
	mat<-read.table(str, header=head,sep=",")
	# Name the different columns
	names(mat)[1]="X"
	names(mat)[2]="Y"
	# If weights weren't specified, fill the weights assuming equal weights for all
	# observations
	if ((is.null(mat[1,3])&&is.null(mat[2,3]))||(mat[1,3]==0&&mat[2,3]==0))
		mat<-calcWeights(mat)
	# Name the weights column
	names(mat)[3]="Weights"
	# Sort the matrix by the X values
	mat<-mat[order(mat[1]),]
	return (mat)
}

# This function finds all the equal X values and replaces them with X value with a Y value
# that is the weighted average of the previous values and with a weight that is equal to
# the sum of all the weights
averageSameXs<-function(mat){
	n<-size(mat)
	temp<-matrix(NA,nrow=n,ncol=3)
	line<-1
	i<-1
	while (i<=n){
		j<-i
		sumY<-0
		sumWeights<-0
		# While the X value is constant (observations with the same X value) sum the
		# observations' Y value and weights
		while (j<=n&mat[i,1]==mat[j,1]){
			sumY<-sumY+mat[j,2]*mat[j,3]
			sumWeights<-sumWeights+mat[j,3]
			j<-j+1
		}
		# Into a temporary matrix fill the Y value and the weight for all the
		# observations that were found with the same X value
		temp[line,1]<-mat[i,1]
		temp[line,2]<-sumY/sumWeights
		temp[line,3]<-sumWeights
		line<-line+1
		i<-j
	}
	# Convert the matrix to a data frame
	temp<-as.data.frame(temp[1:(line-1),])
	# Name the columns accordingly
	names(temp)[1]="X"
	names(temp)[2]="Y"
	names(temp)[3]="Weights"
	return (temp)
}

# This function receives a data set and plot an ACC and LMA graphs for it
# The data source can be manual entry of the data or a csv file (manual entry is the defualt
# if no path is provided).
plotLMA<-function(str=NULL,header=FALSE){
	# If a path to a CSV file wasn't specified, open the manual entry window by calling
	# the fillData function, else call the fillCSVData function
	if (is.null(str)){
		mat<-fillData()
	}
	else{
		mat<-fillCSVData(str,header)
	}
	mat<-calcWeights(mat)
	mat<-averageSameXs(mat)
	mat<-calcFX(mat)
	mat<-calcFY(mat)
	mat<-calcLOI(mat)
	mat<-calcLMA(mat)
	plotGraphs(mat)
	return (mat)
}