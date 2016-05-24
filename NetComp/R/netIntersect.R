#################################
##Name: netIntersect
##Description: Function to find the intersection of 2 square matricies
##O/S: for R
##Date: 10/26/2010 (updated 3/2/2012 to impove memory usage and speed)
##Author: Shannon M. Bell
##Company: Michigan State University
###################################
#This code will take as input (required) 2 square matrixies that must have row/column names
#it can also accept a cutoff (default=NULL) which is the cut off for determining what edges to keep
#note as of 12/13/10 it is done pre-intersection!
#the output will be the average value of edges from input matrix
#Note that the absolute value of the matrix is being used in the calculations.
#if you do not want the absolute value you must change it- and this has potential of removing
#edges that may be highly correlated, just negativly correlated
#on 12/13/2010 altered code to ensure that the edges returned come from both graphs
#change on 1/6/2011 to fix issues with not setting NA to 0 and to try and speed up the process
#3/2/12 using matrix_threshold function to speed up the process

netIntersect<-function(matrix1, matrix2, cutoff=NULL, absolute=TRUE, ...){
    if((row.names(matrix1) != colnames(matrix1)) || (nrow(matrix1) != ncol(matrix1))){
        print("Matrix1 must be square and have same row/column names")
        return(0)
    }
    if((row.names(matrix2) != colnames(matrix2)) || (nrow(matrix2) != ncol(matrix2))){
        print("Matrix2 must be square and have same row/column names")
        return(0)
    }
    #need to put in error for in matrix1 or 2 is not square or is null
    #ordering the graph
    matrix1<-matrix1[sort(row.names(matrix1)), sort(colnames(matrix1))]
    matrix2<-matrix2[sort(row.names(matrix2)), sort(colnames(matrix2))]
    shared.names<-intersect(colnames(matrix1), colnames(matrix2))
    matrix1b<-matrix1[shared.names, shared.names]
    matrix2b<-matrix2[shared.names, shared.names]
    if(is.null(cutoff)){
        cutoff<-0
    }
    #this removes 'data' from lines below cutoff.
    #if cutoff is null then it is 0 and nothing is removed
    matrix1b<-matrix_threshold(matrix1b, threshold=cutoff, minval=0, abs=TRUE, rm.na=TRUE, maxval=NULL)
    matrix2b<-matrix_threshold(matrix2b, threshold=cutoff, minval=0, abs=TRUE, rm.na=TRUE, maxval=NULL)
    
    if(absolute==TRUE){
        matrix3<-(abs(matrix1b) + abs(matrix2b))/2
        matrix3<-matrix_threshold(matrix3, threshold=cutoff, minval=0, abs=TRUE, rm.na=TRUE, maxval=NULL)
    }
    else{
        matrix3<-(matrix1b + matrix2b)/2
        matrix3<-matrix_threshold(matrix3, threshold=cutoff, minval=0, abs=TRUE, rm.na=TRUE, maxval=NULL)
    }
    matrix3
}
