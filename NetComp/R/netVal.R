###########################################################################
# NetVal.c
#
# Authors: 	Shannon M. Bell
#          	Michigan State University
#          	Department of Biochemistry & Molecular Biology
#		Quantitative Biology
#
# Version 1.2	March 6, 2012
#
# Notes: SMB conceived of and wrote the algorithm.
# This algorithm is written to compare 2 adjacency matricies (correlation) and calculate
# the number of TP, FP, TN, FN based on communties defined from a heirctical tree
#matrix1 is the true adj matrix, where matrix2 is the test
#method is the heirictical clustering method
#and h=height, k=#communities passed to cutree
#modified to include the adjusted rand index and balanced accuracy
#################################################3

netVal<-function(matrix1,matrix2, method='ward', k=200, h=NULL, ...){
    if(nrow(matrix1)==ncol(matrix1) && nrow(matrix2)==ncol(matrix2)){
        tree1<-cutree(hclust(as.dist(1-matrix1), method=method), k, h)
        tree2<-cutree(hclust(as.dist(1-matrix2), method=method), k, h)
        #send the trees into netVal to convert the clusters into adjacency matrix
        #the output is a binary graph, where 1 means they cluster together
        x<-.Call("netVal", as.numeric(tree1), pkg="NetComp")
        y<-.Call("netVal",as.numeric(tree2), pkg="NetComp")
        
        row.names(x)<-colnames(matrix1); colnames(x)<-colnames(matrix1)
        row.names(y)<-colnames(matrix2); colnames(y)<-colnames(matrix2)
        #for the next step i need them to have different values
        y<-3*y
        shared.names<-intersect(colnames(x), colnames(y))
        x2<-x[shared.names, shared.names]
        y2<-y[shared.names, shared.names]
        z<-x2+y2
        lz<-lowerTriangle(z, diag=FALSE)
        tp<-length(which(lz==4))
        tn<-length(which(lz==0))
        fp<-length(which(lz==3))
        fn<-length(which(lz==1))
        ##this comuptes the rand index AKA accuracy
        ##the data1 is views as the true or complete dataset
        ##whereas data2 is the imputed or test set
        #RI<-(tp+tn)/(tp+tn+fp+fn)
        #output<-c(TruePos=tp, TrueNeg=tn, FalsePos=fp, FalseNeg=fn, Acc=RI)
        #adding balanced accuracy
        BA<-((.5*tp)/(tp+fn)) + ((.5*tn)/(tn+fp))
        M<-tp+tn+fp+fn
        m1<-tp+fp
        m2<-tp+fn
        EI<-(1/M)*m1*m2
        ARI<-(tp-EI)/(((m1+m2)/2)-EI)
        output<-c(TruePos=tp, TrueNeg=tn, FalsePos=fp, FalseNeg=fn, ARI=ARI, BA=BA)
        output
    }
    else{
        warning("Input matrices must be square")
    }
}