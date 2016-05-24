#################################
##Name: netUnion
##Description: Function to find the union of 2 square matricies (comming out of correlation)
##O/S: for R
##Date: 12/9/2010
##Author: Shannon M. Bell
##Company: Michigan State University
###################################
#This code will take as input (required) 2 square matrixies that must have row/column names
#if you are interested in doing a subset, take those nodes from the corr matrix (still needs to be square)
#the output is a graph with nodes from both party,
#and where a 0.5 is returned if the edge is in 1 graph
# a 1 if in both graphs
#0 if in neither graph (typically due to NA or not meeting cuttoff)
#cutoff represents the cutoff of correlation that should be considered for the intersection
#the abs(input) is used

netUnion<-function(matrix1, matrix2, cutoff=NULL, ...){
    if((row.names(matrix1) != colnames(matrix1)) || (nrow(matrix1) != ncol(matrix1))){
        print("Matrix1 must be square and have same row/column names")
        return(0)
    }
    if((row.names(matrix2) != colnames(matrix2)) || (nrow(matrix2) != ncol(matrix2))){
        print("Matrix2 must be square and have same row/column names")
        return(0)
    }
    #ordering the graph
    matrix1<-matrix1[sort(row.names(matrix1)), sort(colnames(matrix1))]
    matrix2<-matrix2[sort(row.names(matrix2)), sort(colnames(matrix2))]
    #this part transforms mat1 and mat2 into a 1/0 graph
    if(is.null(cutoff)){
        matrix1<-matrix_threshold(matrix1, threshold=0, minval=0, abs=TRUE, rm.na=TRUE, maxval=1)
        matrix2<-matrix_threshold(matrix2, threshold=0, minval=0, abs=TRUE, rm.na=TRUE, maxval=1)
    }
    if(!is.null(cutoff)){
        matrix1<-matrix_threshold(matrix1, threshold=cutoff, minval=0, abs=TRUE, rm.na=TRUE, maxval=1)
        matrix2<-matrix_threshold(matrix2, threshold=cutoff, minval=0, abs=TRUE, rm.na=TRUE, maxval=1)
    }
    
    #need to put in error for in matrix1 or 2 is not square or is null
    shared.names<-intersect(colnames(matrix1), colnames(matrix2))
    #this is the total names, and the length of the final product
    total.names<-sort(union(colnames(matrix1), colnames(matrix2)))
    #need to get the missing from each
    missing2.names<-setdiff(colnames(matrix1), colnames(matrix2))
    missing1.names<-setdiff(colnames(matrix2), colnames(matrix1))
    #get the number of samples
    l1<-length(colnames(matrix1))
    l2<-length(colnames(matrix2))
    #this goes through and addes rows/col of 0 for missing in matrix2
    if(length(missing2.names) >0){
        temp <- matrix(0, nrow = l2, ncol = length(missing2.names))
        temp<-as.data.frame(temp)
        colnames(temp)<-missing2.names
        m2.temp<-cbind(matrix2, temp)
        temp2 <- matrix(0, nrow = length(missing2.names), ncol = length(total.names))
        temp2<-as.data.frame(temp2)
        colnames(temp2)<-colnames(m2.temp)
        rownames(temp2)<-missing2.names
        m2<-rbind(m2.temp, temp2)
    } else{
        m2<-matrix2
    }
    #this goes through and addes rows/col of 0 for missing in matrix1
    if(length(missing1.names) >0){
        temp <- matrix(0, nrow = l1, ncol = length(missing1.names))
        temp<-as.data.frame(temp)
        colnames(temp)<-missing1.names
        m1.temp<-cbind(matrix1, temp)
        temp2 <- matrix(0, nrow = length(missing1.names), ncol = length(total.names))
	temp2 <- as.data.frame(temp2)
	colnames(temp2) <- colnames(m1.temp)
        rownames(temp2)<-missing1.names
        m1<-rbind(m1.temp, temp2)
    } else{
        m1<-matrix1
    }
    matrix1b<-m1[total.names, total.names]
    matrix2b<-m2[total.names, total.names]
    #now i want to combine the 2 graphs
    #result will be a 1/2 if only in 1 matrix and a 1 if in both
    matrix3<-(matrix1b + matrix2b)/2
    matrix3
}
