
#####################################################################################
#                      Functions for Bilder & Loughin (2004)                        #
#####################################################################################

.onLoad<-function(libname, pkgname) {
  assign("MRCV_globals", new.env(), envir=parent.env(environment()))
}

# item.response.table() 
# A function that converts a raw data set to an item response table 
# or data frame (create.dataframe = TRUE)

item.response.table<-function(data, I, J, K = NULL, create.dataframe = FALSE) {
  op<-options()
  on.exit(options(op))
  options(warn=1)
  if((class(data)!="data.frame")&(class(data)!="matrix")) {
    stop("The \"data\" argument requires an object of class \"data.frame\".")
  }
  data<-as.data.frame(data)
  if ((I==0)|(J==0)) {
    stop("\"I\" and \"J\" must be greater than 0.")
  }
  if (!is.numeric(I)|!is.numeric(J)) {
    stop("\"I\" and \"J\" must be numeric.")
  }
  if ((I%%1!=0)|(J%%1!=0)) {
    stop("\"I\" and \"J\" must be integers.")
  }
  I<-as.integer(I)
  J<-as.integer(J)
  if (!is.null(K)) {
    if (!is.numeric(K)) {
      stop("\"K\" must be numeric or NULL.")
    }
    if (K%%1!=0) {
      stop("\"K\" must be an integer.")
    }
    K<-as.integer(K)
    if (K==0) {
      K<-NULL
    }
  }
  if ((class(create.dataframe)!="logical")&(create.dataframe!=1)&(create.dataframe!=0)) {
    warning("The \"create.dataframe\" argument requires a logical value. \n  The input value has been changed to the default value of FALSE.")
    create.dataframe<-FALSE
  }
  create.dataframe<-as.logical(create.dataframe)
  nvars<-1+((I>1)&(J>1))+is.numeric(K)
  if (nvars==1) {
    srcv<-ifelse(I==1,1,I+1)
    mrcv<-ifelse(I==1,2,1)
    r<-length(levels(as.factor(data[,srcv])))
    c<-ifelse(I==1,J,I)
    wy.count<-as.data.frame(matrix(data = NA, nrow = 2*r*c, ncol = 4))
    counter<-0
    for (i in 1:c) {
      wy.count[((counter*r+1):(counter*r+r)),1]<-levels(as.factor(data[,srcv]))
      wy.count[((counter*r+1):(counter*r+r)),2]<-rep(names(data)[(mrcv+i-1)],r)
      wy.count[((counter*r+1):(counter*r+r)),3]<-rep(0,r)
      wy.count[((counter*r+1):(counter*r+r)),4]<-(as.matrix(table(data[,srcv],
                                                 data[,(mrcv+i-1)]))[,1])
      wy.count[((counter*r+r+1):(counter*r+2*r)),1]<-levels(as.factor(data[,srcv]))
      wy.count[((counter*r+r+1):(counter*r+2*r)),2]<-rep(names(data)[(mrcv+i-1)],r)
      wy.count[((counter*r+r+1):(counter*r+2*r)),3]<-rep(1,r)
      wy.count[((counter*r+r+1):(counter*r+2*r)),4]<-(as.matrix(table(data[,srcv],
                                                     data[,(mrcv+i-1)]))[,2])
      counter<-counter+2
    }
    wy.count[,1]<-factor(wy.count[,1])
    wy.count[,2]<-factor(wy.count[,2],names(data)[mrcv:(mrcv+c-1)])
    wy.count[,3]<-as.numeric(wy.count[,3])
    wy.count[,4]<-as.numeric(wy.count[,4])
    colnames(wy.count) <- c("W", "Y", "yj", "count")
    if (!create.dataframe) {
      output<-tabular(Heading()*count*Heading()*W
                        ~Heading()*Y*Heading()*Factor(yj, yj, c("0","1   "))*Heading()*mean, 
                        data=wy.count)
    }
    if (create.dataframe) {
      output<-wy.count
    }
  }
  nrows<-(2^nvars)*I*J*max(1,K)
  if (nvars==2) {
    wy.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (2*nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        counter<-counter+2^nvars
        cell.count<-table(data[,i],data[,I+j])
        wy.count[(counter-3),]<-c(names(data)[i], names(data)[(I+j)], 0, 0, 
                                  cell.count[1,1])
        wy.count[(counter-2),]<-c(names(data)[i], names(data)[(I+j)], 0, 1, 
                                  cell.count[1,2])
        wy.count[(counter-1),]<-c(names(data)[i], names(data)[(I+j)], 1, 0, 
                                  cell.count[2,1])
        wy.count[(counter),]  <-c(names(data)[i], names(data)[(I+j)], 1, 1, 
                                  cell.count[2,2])
      }
    }
    wy.count[,1]<-factor(wy.count[,1],names(data)[1:I])
    wy.count[,2]<-factor(wy.count[,2],names(data)[(I+1):(I+J)])
    wy.count[,(nvars+1):(2*nvars+1)]<-lapply(wy.count[,(nvars+1):(2*nvars+1)], 
                                             as.numeric)
    colnames(wy.count) <- c("W", "Y", "wi", "yj", "count")
    if (!create.dataframe) {
      output<-tabular(count*(mean)*W*Heading()*Factor(wi, wi, c("0","1 "))
                        ~Heading()*Y*Heading()*Factor(yj, yj, c("0","1   ")), data=wy.count, 
                        suppressLabels=3)
    }
    if (create.dataframe) {
      output<-wy.count
    }
  }
  if (nvars==3) {
    wyz.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (2*nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          counter<-counter+2^nvars
          cell.count<-table(data[,i],data[,(I+j)],data[,(I+J+k)])
          wyz.count[(counter-7),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,0,0,cell.count[1,1,1])
          wyz.count[(counter-6),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,0,1,cell.count[1,1,2])
          wyz.count[(counter-5),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,1,0,cell.count[1,2,1])
          wyz.count[(counter-4),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],0,1,1,cell.count[1,2,2])
          wyz.count[(counter-3),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,0,0,cell.count[2,1,1])
          wyz.count[(counter-2),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,0,1,cell.count[2,1,2])
          wyz.count[(counter-1),]<-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,1,0,cell.count[2,2,1])
          wyz.count[(counter),]  <-c(names(data)[i],names(data)[(I+j)],
                                     names(data)[(I+J+k)],1,1,1,cell.count[2,2,2])
        }
      }
    }
    wyz.count[,1]<-factor(wyz.count[,1],names(data)[1:I])
    wyz.count[,2]<-factor(wyz.count[,2],names(data)[(I+1):(I+J)])
    wyz.count[,3]<-factor(wyz.count[,3],names(data)[(I+J+1):(I+J+K)])
    wyz.count[,(nvars+1):(2*nvars+1)]<-lapply(wyz.count[,(nvars+1):(2*nvars+1)], 
                                              as.numeric)
    colnames(wyz.count) <- c("W", "Y", "Z", "wi", "yj", "zk", "count")
    if (!create.dataframe) {
      output<-vector("list", K*2)
      output.names<-matrix(data=NA, nrow=1, ncol=K*2)
      counter<-0
      for (k in 1:K) {
        for (c in 0:1) {
          counter<-counter+1
          output[counter]<-list(tabular(count*(mean)*W*Heading()*Factor(wi, wi, c("0","1 "))
                            ~Heading()*Y*Heading()*Factor(yj, yj, c("0","1   ")), 
                            data=wyz.count[((wyz.count[,3]==names(data)[I+J+k])&(wyz.count[,6]==c)),], 
                            suppressLabels=3))
          output.names[1,counter]<-paste(names(data[(I+J+k)]), "=", c)
        }
      }
      names(output)<-output.names
    }
    if (create.dataframe) {
      output<-wyz.count
    }
  }
  output
}

#####################################################################################

# marginal.table() 
# A function that converts a raw data set to a marginal table

marginal.table<-function(data, I, J, K = NULL) {
  op<-options()
  on.exit(options(op))
  options(warn=1)
  if((class(data)!="data.frame")&(class(data)!="matrix")) {
    stop("The \"data\" argument requires an object of class \"data.frame\".")
  }
  data<-as.data.frame(data)
  if ((I==0)|(J==0)) {
    stop("\"I\" and \"J\" must be greater than 0.")
  }
  if (!is.numeric(I)|!is.numeric(J)) {
    stop("\"I\" and \"J\" must be numeric.")
  }
  if ((I%%1!=0)|(J%%1!=0)) {
    stop("\"I\" and \"J\" must be integers.")
  }
  I<-as.integer(I)
  J<-as.integer(J)
  if (!is.null(K)) {
    if (!is.numeric(K)) {
      stop("\"K\" must be numeric or NULL.")
    }
    if (K%%1!=0) {
      stop("\"K\" must be an integer.")
    }
    K<-as.integer(K)
    if (K==0) {
      K<-NULL
    }
  }
  nvars<-1+((I>1)&(J>1))+is.numeric(K)
  if (nvars==1) {
    srcv<-ifelse(I==1,1,I+1)
    mrcv<-ifelse(I==1,2,1)
    r<-length(levels(as.factor(data[,srcv])))
    c<-ifelse(I==1,J,I)
    wy.count<-as.data.frame(matrix(data = NA, nrow = r, ncol = 4))
    merge.matrix<-as.data.frame(matrix(data= NA, nrow = 1, ncol = 4))
    counter<-0
    for (i in 1:c) {
      wy.count[1:r,1]<-levels(as.factor(data[,srcv]))
      wy.count[1:r,2]<-rep(names(data)[(mrcv+i-1)],r)
      wy.count[1:r,3]<-as.matrix(table(data[,srcv],data[,(mrcv+i-1)]))[,2]
      wy.count[1:r,4]<-as.matrix(table(data[,srcv]))
      if (i==1) {
        for (j in 1:r) {
          rep.matrix<-matrix(as.vector(rep(as.vector(wy.count[j,]),
                             wy.count[j,4]-(c-1))), (wy.count[j,4]-(c-1)), 4, 
                             byrow=TRUE)
          merge.matrix<-rbind(merge.matrix, rep.matrix)
        }
        n.iplus<-as.matrix(table(data[,srcv]))
      }
      if (i > 1) {
        merge.matrix<-rbind(merge.matrix,wy.count)
      }
      counter<-counter+1
    }
    merge.matrix<-merge.matrix[2:nrow(merge.matrix),]
    merge.matrix[,1]<-factor(unlist(merge.matrix[,1]))
    merge.matrix[,2]<-factor(unlist(merge.matrix[,2]),names(data)[mrcv:(mrcv+c-1)])
    merge.matrix[,3]<-as.numeric(unlist(merge.matrix[,3]))
    merge.matrix[,4]<-as.numeric(unlist(merge.matrix[,4]))
    colnames(merge.matrix) <- c("W", "Y", "count", "n")
    MRCV_globals$n.iplus<-n.iplus
    MRCV_globals$c<-c
    MRCV_globals$counter<--1
    calc.percent <- function(x) {
      MRCV_globals$counter<-MRCV_globals$counter+1
      perc<-round(100*mean(x)/MRCV_globals$n.iplus[(trunc(MRCV_globals$counter/MRCV_globals$c)+1)], 2)
      if (perc==0) {
        perc<-"0.00"
      }
      perc
    }
    marginal.table<-tabular(Heading()*count*Heading()*W
                            ~Justify(l)*(n=1)+Heading()*Y*(Heading("count")*mean + 
                            Heading("  %     ")*calc.percent), data=merge.matrix)
  }
  n<-nrow(data)
  nrows<-I*J*max(1,K)
  calc.percent <- function(x) {
    perc<-round(100*mean(x)/n, 2)
    if (perc==0) {
      perc<-"0.00"
    }
    perc
  }
  if (nvars==2) {
    wy.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        counter<-counter+1
        cell.count<-table(data[,i],data[,(I+j)])
        wy.count[counter,]<-c(names(data)[i],names(data)[(I+j)],cell.count[2,2]) 
      }
    }
    wy.count[,1]<-factor(wy.count[,1],names(data)[1:I])
    wy.count[,2]<-factor(wy.count[,2],names(data)[(I+1):(I+J)])
    wy.count[,(nvars+1)]<-as.numeric(wy.count[,(nvars+1)])
    colnames(wy.count) <- c("W", "Y", "count")
    marginal.table<-tabular(count*W~Heading()*Y*(Heading("count")*mean + 
                            Heading("  %     ")*calc.percent), data=wy.count, 
                            suppressLabels=2)
  }
  if (nvars==3) {
    wyz.count<-as.data.frame(matrix(data = NA, nrow = nrows, ncol = (nvars+1)))
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
        counter<-counter+1
        cell.count<-table(data[,i],data[,(I+j)],data[,(I+J+k)])
        wyz.count[counter,]<-c(names(data)[i],names(data)[(I+j)],names(data)[(I+J+k)],
                               cell.count[2,2,2]) 
        }
      }
    }
    wyz.count[,1]<-factor(wyz.count[,1],names(data)[1:I])
    wyz.count[,2]<-factor(wyz.count[,2],names(data)[(I+1):(I+J)])
    wyz.count[,3]<-factor(wyz.count[,3],names(data)[(I+J+1):(I+J+K)])
    wyz.count[,(nvars+1)]<-as.numeric(wyz.count[,(nvars+1)])
    colnames(wyz.count) <- c("W", "Y", "Z", "count")
    marginal.table<-vector("list", K)
    for (k in 1:K) {
      marginal.table[k]<-list(tabular(count*W~Heading()*Y*(Heading("count")*mean + 
                              Heading("  %     ")*calc.percent), 
                              data=wyz.count[(wyz.count[,3]==names(data)[I+J+k]),], 
                              suppressLabels=2))
    }
    names(marginal.table)<-names(data[(I+J+1):(I+J+K)])
  }
  marginal.table
}

#####################################################################################

# check.zero()
# Function that adds a small constant to 0 cell counts

check.zero<-function(value, add.constant = .5) {
  if (value==0) newvalue<-add.constant else newvalue<-value
  newvalue
}

#####################################################################################

# MMI.stat()
# A function that calculates X^2_S and X^2_S.ij to test for MMI (1 MRCV case)
# Called by the general MI.stat() function

MMI.stat<-function(data, I, J, summary.data, add.constant) {
  op<-options()
  on.exit(options(op))
  srcv<-ifelse(I==1,1,I+1)
  mrcv<-ifelse(I==1,2,1)
  c<-ifelse(I==1,J,I)
  #For inputted data set that is a summary file
  if (summary.data) {
    r<-length(levels(as.factor(data[,1])))
    data[,1:2]<-lapply(data[,1:2], factor)
    X.sq.S.ij<-matrix(data = NA, nrow = 1, ncol = c)
    counter<-0
    for (i in 1:c) {
      table.11<-data[((data[,2]==levels(data[,2])[i])&(data[,3]==0)), 4]
      table.12<-data[((data[,2]==levels(data[,2])[i])&(data[,3]==1)), 4]
      n.table<-matrix(data = c(table.11, table.12), nrow = r, ncol = 2)
      #Only calculate statistics for valid rx2 tables
      if (min(c(rowSums(n.table), colSums(n.table))) > 0) {
        counter<-counter+1
        #Add .5 to 0 cell counts
        n.table<-apply(X = n.table, MARGIN = c(1,2), FUN = check.zero, 
                       add.constant = add.constant)
        options(warn = -1) 
        X.sq.S.ij[1,i]<-chisq.test(n.table,correct=F)$statistic
        options(warn = 0) 
      }
    }
    colnames(X.sq.S.ij)<-levels(data[,2])
    output<-list(X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins = 
                   counter)
  }  
  
  #For inputted data set that is a raw file
  if (!summary.data) {
    r<-length(levels(as.factor(data[,srcv])))
    X.sq.S.ij<-matrix(data = NA, nrow = 1, ncol = c)
    counter<-0
    for (i in 1:c) {
      #Only calculate statistics for valid tables (correct dim = rx2)
      if (sum(dim(table(data[,srcv],data[,(mrcv+i-1)])))==(r+2)){
        counter<-counter+1
        n.table<-table(data[,srcv],data[,(mrcv+i-1)])
        #Add .5 to 0 cell counts
        n.table<-apply(X = n.table, MARGIN = c(1,2), FUN = check.zero, 
                       add.constant = add.constant)
        options(warn = -1) 
        X.sq.S.ij[1, i]<-chisq.test(n.table,correct=F)$statistic
        options(warn = 0) 
      }
    }
    colnames(X.sq.S.ij)<-names(data)[mrcv:(mrcv+c-1)]
    output<-list(X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins = 
                   counter)
  }
  output
}

#####################################################################################

# SPMI.stat()
# A function that calculates X^2_S and X^2_S.ij to test for SPMI (2 MRCV case)
# Called by the general MI.stat() function

SPMI.stat<-function(data, I, J, summary.data, add.constant) {
  op<-options()
  on.exit(options(op))
  #For inputted data set that is a summary file
  if (summary.data) {
    data[,1:2]<-lapply(data[,1:2], factor)
    X.sq.S.ij<-matrix(data = NA, nrow = I, ncol = J)
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        table.11<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==0)&(data[,4]==0)),5]
        table.12<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==0)&(data[,4]==1)),5]
        table.21<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==1)&(data[,4]==0)),5]
        table.22<-data[((data[,1]==levels(data[,1])[i])& 
                      (data[,2]==levels(data[,2])[j])&(data[,3]==1)&(data[,4]==1)),5]
        n.table<-matrix(data = c(table.11, table.21, table.12, table.22), 
                        nrow = 2, ncol = 2)
        #Only calculate statistics for valid 2x2 tables (no 00 in row or col)
        if (min(c(rowSums(n.table), colSums(n.table))) > 0) {
          counter<-counter+1
          #Add .5 to 0 cell counts
          n.table<-apply(X = n.table, MARGIN = c(1,2), FUN = check.zero, 
                         add.constant = add.constant)
          options(warn = -1) 
          X.sq.S.ij[i,j]<-chisq.test(n.table,correct=F)$statistic
          options(warn = 0) 
        }
      }
    }
    rownames(X.sq.S.ij)<-levels(data[,1])
    colnames(X.sq.S.ij)<-levels(data[,2])
    output<-list(X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins = 
                   counter)
  }  
  
  #For inputted data set that is a raw file
  if (!summary.data) {
    X.sq.S.ij<-matrix(data = NA, nrow = I, ncol = J)
    counter<-0
    for (i in 1:I) {
      for (j in 1:J) {
        #Only calculate statistics for valid tables (correct dim = 2x2)
        if (sum(dim(table(data[,i],data[,(I+j)])))==4){
          counter<-counter+1
          n.table<-table(data[,i],data[,(I+j)])
          #Add .5 to 0 cell counts
          n.table<-apply(X = n.table, MARGIN = c(1,2), FUN = check.zero, 
                         add.constant = add.constant)
          options(warn = -1)
          X.sq.S.ij[i, j]<-chisq.test(n.table,correct=F)$statistic
          options(warn = 0) 
        }
      }
    }
    rownames(X.sq.S.ij)<-names(data)[1:I]
    colnames(X.sq.S.ij)<-names(data)[(I+1):(I+J)]
    output<-list(X.sq.S = sum(X.sq.S.ij), X.sq.S.ij = X.sq.S.ij, valid.margins = 
                   counter)
  }
  output
}

#####################################################################################

# MI.stat()
# A function that calculates X^2_S and X^2_S.ij

MI.stat<-function(data, I, J, summary.data = FALSE, add.constant = .5) {
  nvars<-1+((I>1)&(J>1))
  if (nvars==1) {
    output<-MMI.stat(data = data, I = I, J = J, summary.data = summary.data, 
                     add.constant = add.constant)
  }
  if (nvars==2) {
    output<-SPMI.stat(data = data, I = I, J = J, summary.data = summary.data, 
                     add.constant = add.constant)
  }
  output
}

#####################################################################################

# check.min()
# Function used in Bonferroni adjustments (p-value should not be greater than 1)

check.min<-function(value) {
  min(value, 1.00000)
}

#####################################################################################

# print.MMI()
# Function used to control display of output with class MMI

print.MMI<-function(x, ...) {
  op<-options()
  on.exit(options(op))
  options(scipen=10)
  type<-names(x)
  data<-x$general$data
  I<-x$general$I
  J<-x$general$J
  summary.data<-x$general$summary.data
  X.sq.S<-round(x$general$X.sq.S, 2)
  X.sq.S.ij<-round(x$general$X.sq.S.ij, 2)
  srcv<-ifelse(I==1,1,I+1)
  mrcv<-ifelse(I==1,2,1)
  c<-ifelse(I==1,J,I)
  if (summary.data) {
    r<-length(levels(as.factor(data[,1])))
    rownames(X.sq.S.ij)<-""
    colnames(X.sq.S.ij)<-levels(data[,2]) 
  }
  if (!summary.data) {
    r<-length(levels(as.factor(data[,srcv])))
    rownames(X.sq.S.ij)<-""
    colnames(X.sq.S.ij)<-names(data)[mrcv:(mrcv+c-1)]
  }
  cat("Test for Multiple Marginal Independence (MMI)", "\n", "\n")
  cat("Unadjusted Pearson Chi-Square Tests for Independence:", "\n")
  cat("X^2_S =", X.sq.S, "\n")
  cat("X^2_S.ij =", "\n") 
  print.default(X.sq.S.ij)
  cat("\n")
  if (any(type=="boot")){
    B.discard<-x$boot$B.discard
    B.use<-x$boot$B.use
    p.value.boot<-as.name(paste("=", round(x$boot$p.value.boot, 4))) 
    if (x$boot$p.value.boot == 0) {
      p.value.boot<-as.name(paste("<", round(1/B.use, 4)))
    }
    p.combo.prod<-as.name(paste("=", round(x$boot$p.combo.prod.boot, 4)))
    if (x$boot$p.combo.prod.boot == 0) {
      p.combo.prod<-as.name(paste("<", round(1/B.use, 4)))
    }
    p.combo.min<-as.name(paste("=", round(x$boot$p.combo.min.boot, 4)))
    if (x$boot$p.combo.min.boot ==0) {
      p.combo.min<-as.name(paste("<", round(1/B.use, 4)))
    }
    cat("Bootstrap Results:", "\n")
    if (B.discard > 0) {
      cat(B.discard, "resamples were removed from the analysis due to")  
      cat(" not having all rows or columns represented in an rx2 table", "\n")
    }
    cat("Final results based on", B.use, "resamples", "\n")
    cat("p.boot", p.value.boot, "\n")
    cat("p.combo.prod", p.combo.prod, "\n")
    cat("p.combo.min", p.combo.min, "\n", "\n")
  }
  if (any(type=="rs2")){
    X.sq.S.rs2<-round(x$rs2$X.sq.S.rs2, 2)
    df.rs2<-round(x$rs2$df.rs2, 2)
    p.value.rs2<-as.name(paste("=", round(x$rs2$p.value.rs2, 4)))
    if (round(x$rs2$p.value.rs2, 4) < .0001) {
      p.value.rs2<-as.name(paste("<", .0001))
    }
    cat("Second-Order Rao-Scott Adjusted Results:", "\n")
    cat("X^2_S.adj =", X.sq.S.rs2, "\n")
    cat("df.adj =", df.rs2, "\n")
    cat("p.adj", p.value.rs2, "\n", "\n")
  }
  if (any(type=="bon")){
    p.value.bon<-as.name(paste("=", round(x$bon$p.value.bon, 4)))
    if (round(x$bon$p.value.bon, 4) < .0001) {
      p.value.bon<-as.name(paste("<", .0001))
    }
    X.sq.S.ij.p.bon<-format(round(x$bon$X.sq.S.ij.p.bon,4), 
                           digits = 4)
    if (summary.data) {
      rownames(X.sq.S.ij.p.bon)<-""
      colnames(X.sq.S.ij.p.bon)<-levels(data[,2]) 
    }
    if (!summary.data) {
      rownames(X.sq.S.ij.p.bon)<-""
      colnames(X.sq.S.ij.p.bon)<-names(data)[mrcv:(mrcv+c-1)] 
    }
    cat("Bonferroni Adjusted Results:", "\n")
    cat("p.adj", p.value.bon, "\n")    
    cat("p.ij.adj =", "\n")
    print.default(X.sq.S.ij.p.bon, quote=FALSE)
    cat("\n")
  } 
  invisible(x)
}

#####################################################################################

# print.SPMI()
# Function used to control display of output with class SPMI

print.SPMI<-function(x, ...) {
  op<-options()
  on.exit(options(op))
  options(scipen=10)
  type<-names(x)
  data<-x$general$data
  I<-x$general$I
  J<-x$general$J
  summary.data<-x$general$summary.data
  X.sq.S<-round(x$general$X.sq.S, 2)
  X.sq.S.ij<-round(x$general$X.sq.S.ij, 2)
  if (summary.data) {
    rownames(X.sq.S.ij)<-levels(data[,1])
    colnames(X.sq.S.ij)<-levels(data[,2]) 
  }
  if (!summary.data) {
    rownames(X.sq.S.ij)<-names(data)[1:I]
    colnames(X.sq.S.ij)<-names(data)[(I+1):(I+J)]
  }
  cat("Test for Simultaneous Pairwise Marginal Independence (SPMI)", "\n", "\n")
  cat("Unadjusted Pearson Chi-Square Tests for Independence:", "\n")
  cat("X^2_S =", X.sq.S, "\n")
  cat("X^2_S.ij =", "\n") 
  print.default(X.sq.S.ij)
  cat("\n")
  if (any(type=="boot")){
    B.discard<-x$boot$B.discard
    B.use<-x$boot$B.use
    p.value.boot<-as.name(paste("=", round(x$boot$p.value.boot, 4))) 
    if (x$boot$p.value.boot == 0) {
      p.value.boot<-as.name(paste("<", round(1/B.use, 4)))
    }
    p.combo.prod<-as.name(paste("=", round(x$boot$p.combo.prod.boot, 4)))
    if (x$boot$p.combo.prod.boot == 0) {
      p.combo.prod<-as.name(paste("<", round(1/B.use, 4)))
    }
    p.combo.min<-as.name(paste("=", round(x$boot$p.combo.min.boot, 4)))
    if (x$boot$p.combo.min.boot == 0) {
      p.combo.min<-as.name(paste("<", round(1/B.use, 4)))
    }
    cat("Bootstrap Results:", "\n")
    if (B.discard > 0) {
      cat(B.discard, "resamples were removed from the analysis due to")  
      cat(" not having all rows or columns represented in a 2x2 table", "\n")
    }
    cat("Final results based on", B.use, "resamples", "\n")
    cat("p.boot", p.value.boot, "\n")
    cat("p.combo.prod", p.combo.prod, "\n")
    cat("p.combo.min", p.combo.min, "\n", "\n")
  }
  if (any(type=="rs2")){
    X.sq.S.rs2<-round(x$rs2$X.sq.S.rs2, 2)
    df.rs2<-round(x$rs2$df.rs2, 2)
    p.value.rs2<-as.name(paste("=", round(x$rs2$p.value.rs2, 4)))
    if (round(x$rs2$p.value.rs2, 4) < .0001) {
      p.value.rs2<-as.name(paste("<", .0001))
    }
    cat("Second-Order Rao-Scott Adjusted Results:", "\n")
    cat("X^2_S.adj =", X.sq.S.rs2, "\n")
    cat("df.adj =", df.rs2, "\n")
    cat("p.adj", p.value.rs2, "\n", "\n")
  }
  if (any(type=="bon")){
    p.value.bon<-as.name(paste("=", round(x$bon$p.value.bon, 4)))
    if (round(x$bon$p.value.bon, 4) < .0001) {
      p.value.bon<-as.name(paste("<", .0001))
    }
    X.sq.S.ij.p.bon<-format(round(x$bon$X.sq.S.ij.p.bon,4), 
                            digits = 4)
    rownames(X.sq.S.ij.p.bon)<-names(data)[1:I]
    colnames(X.sq.S.ij.p.bon)<-names(data)[(I+1):(I+J)] 
    if (summary.data) {
      rownames(X.sq.S.ij.p.bon)<-levels(data[,1])
      colnames(X.sq.S.ij.p.bon)<-levels(data[,2]) 
    }
    cat("Bonferroni Adjusted Results:", "\n")
    cat("p.adj", p.value.bon, "\n")    
    cat("p.ij.adj =", "\n")
    print.default(X.sq.S.ij.p.bon, quote=FALSE)
    cat("\n")
  }
  invisible(x)
}

#####################################################################################

# MMI.test()
# A function that performs bootstrap testing, the second-order Rao-Scott approach, 
#   and Bonferroni adjustments to test for MMI (1 MRCV case)
# Called by the general MI.test() function

MMI.test<-function(data, I, J, type, B, B.max, summary.data, add.constant, plot.hist, 
                    print.status) {
  n<-nrow(data)
  srcv<-ifelse(I==1,1,I+1)
  mrcv<-ifelse(I==1,2,1)
  c<-ifelse(I==1,J,I)
  if (summary.data) {
    r<-length(levels(as.factor(data[,1])))
  }
  if (!summary.data) {
    r<-length(levels(as.factor(data[,srcv])))
  }
  #Observed statistics 
  observed<-MI.stat(data = data, I = I, J = J, 
                      summary.data = summary.data, add.constant = add.constant)
  p.value.ij<-1 - pchisq(q = observed$X.sq.S.ij, df = (r-1))
  p.value.min<-min(p.value.ij)
  p.value.prod<-prod(p.value.ij)    
  
  #Bootstrap
  if (type == "boot" | type == "all") {
    X.sq.S.star<-numeric(length(B.max))
    p.value.b.min<-numeric(length(B.max))
    p.value.b.prod<-numeric(length(B.max))
    X.sq.S.ij.star<-matrix(data = NA, nrow = c, ncol = B.max)
    discard<-0 #Count number of resamples with incorrect dimensions
    counter<-0
    b<-0
    filled<-0
    # create progress bar
    if (print.status) {
      cat("Bootstrap Progress:", "\n")
      pb <- txtProgressBar(min = 0, max = B.max, style = 3)
    }
    while(((counter>=B)+(b>=B.max)) < 1) {
      b<-b+1
      n_iplus<-as.matrix(table(data[,srcv])) 
      W<-rep(levels(as.factor(data[,srcv])), times = n_iplus)
      Y<-sample(x = 1:n, size = n, replace = TRUE)
      data.star<-cbind(W, data[Y, mrcv:(mrcv+c-1)])
      stat.star<-MI.stat(data = data.star, I = 1, J = c, 
                         summary.data = FALSE, add.constant = add.constant)
      discard<-discard+(stat.star$valid.margins<c) #Discard if any table < rx2
      drop<-1
      if (stat.star$valid.margins==c){ 
        counter<-counter+1
        X.sq.S.star[counter]<-stat.star$X.sq.S
        X.sq.S.ij.star[,counter]<-stat.star$X.sq.S.ij
        p.value.ij<-1 - pchisq(q = stat.star$X.sq.S.ij, df = (r-1))
        p.value.b.min[counter]<-min(p.value.ij)
        p.value.b.prod[counter]<-prod(p.value.ij)
        drop<-0
      }
      # update progress bar
      if (B==B.max) {
        if (print.status) {
          setTxtProgressBar(pb, b)
        }
      }
      if (B!=B.max) {
        if (print.status) {
          left<-B.max-filled
          expect<-(1-drop)*min(left,(B-counter))+drop*max(left,(B.max-b))
          filled<-filled+left/expect
          setTxtProgressBar(pb, filled)
        }
      }
    }
    if (print.status) {
      setTxtProgressBar(pb, B.max)
      close(pb)
    }
    B.use<-min(B,(B.max-discard)) #Only use desired number of resamples
    X.sq.S.star<-X.sq.S.star[1:B.use]
    X.sq.S.ij.star<-X.sq.S.ij.star[,1:B.use]
    p.value.b.min<-p.value.b.min[1:B.use]
    p.value.b.prod<-p.value.b.prod[1:B.use]
    p.value.boot<-mean(X.sq.S.star>=observed$X.sq.S)
    p.combo.min<-list(p = p.value.min, p.star = p.value.b.min, 
                      overall.p = mean(p.value.b.min<=p.value.min))
    p.combo.prod<-list(p = p.value.prod, p.star = p.value.b.prod, 
                       overall.p = mean(p.value.b.prod<=p.value.prod))
    #Histograms
    if (plot.hist) {
      layout.m <- matrix(c(1, 0, 1, 3, 2, 3, 2, 0), nrow = 2, ncol = 4)
      layout(layout.m)
      hist(x = X.sq.S.star, main = "Test using sum statistic", xlab = 
             expression(X[S]^{"2*"}), xlim = c(min(X.sq.S.star, observed$X.sq.S), 
                                               max(X.sq.S.star, observed$X.sq.S)))
      abline(v = observed$X.sq.S, col = "darkgreen", lwd = 5)
      hist(x = p.value.b.prod, main = "Test using product of p-values", xlab = 
             expression(tilde(p)[prod]^{"*"}), xlim = c(min(p.value.b.prod, 
             p.value.prod), max(p.value.b.prod, p.value.prod)))
      abline(v = p.value.prod, col = "darkgreen", lwd = 5)
      hist(x = p.value.b.min, main = "Test using minimum of p-values", xlab = 
             expression(tilde(p)[min]^{"*"}), xlim = c(min(p.value.b.min, 
             p.value.min), max(p.value.b.min, p.value.min)))
      abline(v = p.value.min, col = "darkgreen", lwd = 5)
    }
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                      summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                      observed$X.sq.S.ij), boot = list(B.use = B.use, 
                      B.discard = discard, p.value.boot = p.value.boot, 
                      p.combo.min.boot = p.combo.min$overall.p, p.combo.prod.boot = 
                      p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star, X.sq.S.ij.star = 
                      X.sq.S.ij.star, p.combo.min.star = p.combo.min$p.star, 
                      p.combo.prod.star = p.combo.prod$p.star))
    output.boot<-list(B.use = B.use, B.discard = discard, 
                      p.value.boot = p.value.boot, p.combo.min.boot = 
                      p.combo.min$overall.p, p.combo.prod.boot = 
                      p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star, 
                      X.sq.S.ij.star = X.sq.S.ij.star, p.combo.min.star = 
                      p.combo.min$p.star, p.combo.prod.star = p.combo.prod$p.star)
  }
  
  #Second-order Rao-Scott adjustment
  if (type == "rs2" | type == "all") {
    #Appendix A calculations
    n.counts<-as.data.frame(table(data)) 
    cols <- c(srcv,(mrcv:(mrcv+c-1))) 
    n.counts.ord<-n.counts[do.call("order", as.data.frame(n.counts[,cols])), ]    
    Y.counts<-as.data.frame(table(data[,(mrcv:(mrcv+c-1))])) 
    cols <- c(1:c) 
    Y.counts.ord<-Y.counts[do.call("order", as.data.frame(Y.counts[,cols])), ]
    n_iplus<-as.matrix(table(data[,srcv])) 
    tau<-n.counts.ord[,ncol(n.counts.ord)]/rep(n_iplus, each = 2^c)   
    G.tilde<-t(data.matrix(Y.counts.ord[,1:c])-1)  
    I.r<-diag(r)
    G<-kronecker(I.r,G.tilde)
    pi<-G%*%tau
    m<-pi*rep(n_iplus, each = c)
    a.i<-n_iplus/n
    pi.not.j<-(1/n)*kronecker(matrix(rep(1,r),1,r),diag(c))%*%m
    j.r<-matrix(1,r,1)
    pi.not<-kronecker(j.r,pi.not.j)
    I.rc<-diag(r*c)
    I.c<-diag(c)
    J.rr<-matrix(1,r,r)
    A<-diag(as.vector(a.i))
    H<-I.rc-kronecker(J.rr%*%A,I.c)
    D<-kronecker(diag(as.vector(n/n_iplus)),diag(as.vector(pi.not.j*(1-pi.not.j))))
    V<-matrix(data = 0, nrow = r*2^c, ncol = r*2^c)
    for (i in 1:r) {
       V[((i-1)*2^c+1):((i-1)*2^c+2^c),((i-1)*2^c+1):((i-1)*2^c+2^c)]<-((1/a.i[i])
          *(diag(as.vector(tau[((i-1)*2^c+1):((i-1)*2^c+2^c)]))
          -tcrossprod(tau[((i-1)*2^c+1):((i-1)*2^c+2^c)])))
    }
    Di.HGVGH<-diag(1/diag(D),dim(D))%*%H%*%G%*%tcrossprod(tcrossprod(V,G),H)
    Di.HGVGH.eigen<-Re(eigen(Di.HGVGH)$values) 
    sum.Di.HGVGH.eigen.sq<-sum(Di.HGVGH.eigen^2)
    X.sq.S.rs2<-(r-1)*c*observed$X.sq.S/sum.Di.HGVGH.eigen.sq 
    df.rs2<-(r-1)^2*c^2/sum.Di.HGVGH.eigen.sq       
    X.sq.S.p.value.rs2<-1-pchisq(q = X.sq.S.rs2, df = df.rs2) 
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                 summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                 observed$X.sq.S.ij), rs2 = list(X.sq.S.rs2 = X.sq.S.rs2, 
                 df.rs2 = df.rs2, p.value.rs2 = X.sq.S.p.value.rs2))
    output.rs2<-list(X.sq.S.rs2 = X.sq.S.rs2, df.rs2 = df.rs2, 
                 p.value.rs2 = X.sq.S.p.value.rs2)
  }
  
  #Bonferroni
  if (type == "bon" | type == "all") {   
    p.value.bon<-min(p.value.min*c,1) 
    X.sq.S.ij.p.bon<-apply(X = c*(1 - pchisq(q = observed$X.sq.S.ij, 
                           df = (r-1))), MARGIN = c(1,2), FUN = check.min)
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                 summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                 observed$X.sq.S.ij), bon = list(p.value.bon = p.value.bon, 
                 X.sq.S.ij.p.bon = X.sq.S.ij.p.bon))
    output.bon<-list(p.value.bon = p.value.bon, X.sq.S.ij.p.bon = X.sq.S.ij.p.bon) 
  }
  if (type == "all") {
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                 summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                 observed$X.sq.S.ij), boot = output.boot, rs2 = output.rs2, 
                 bon = output.bon)
  }  
  class(output)<-"MMI"
  output
}

#####################################################################################

# SPMI.test()
# A function that performs bootstrap testing, the second-order Rao-Scott approach, 
#   and Bonferroni adjustments to test for SPMI (2 MRCV case)
# Called by the general MI.test() function

SPMI.test<-function(data, I, J, type, B, B.max, summary.data, add.constant, plot.hist, 
                    print.status) {
  n<-nrow(data)
  #Observed statistics 
  observed<-MI.stat(data = data, I = I, J = J, 
                    summary.data = summary.data, add.constant = add.constant)
  p.value.ij<-1 - pchisq(q = observed$X.sq.S.ij, df = 1)
  p.value.min<-min(p.value.ij)
  p.value.prod<-prod(p.value.ij)    
  
  #Bootstrap
  if (type == "boot" | type == "all") {
    X.sq.S.star<-numeric(length(B.max))
    p.value.b.min<-numeric(length(B.max))
    p.value.b.prod<-numeric(length(B.max))
    X.sq.S.ij.star<-matrix(data=NA, nrow = I*J, ncol = B.max)
    discard<-0 #Count number of resamples with incorrect dimensions
    counter<-0
    b<-0
    filled<-0
    # create progress bar
    if (print.status) {
      cat("Bootstrap Progress:", "\n")
      pb <- txtProgressBar(min = 0, max = B.max, style = 3)
    }
    while(((counter>=B)+(b>=B.max)) < 1) {
      b<-b+1
      #Resample W_s and Y_s independently
      W<-sample(x = 1:n, size = n, replace = TRUE) 
      Y<-sample(x = 1:n, size = n, replace = TRUE)
      data.star<-cbind(data[W, 1:I], data[Y, (I+1):(I+J)])
      stat.star<-MI.stat(data = data.star, I = I, J = J, 
                         summary.data = FALSE, add.constant = add.constant)
      discard<-discard+(stat.star$valid.margins<I*J) #Discard if any table < 2x2
      drop<-1
      if (stat.star$valid.margins==I*J){ 
        counter<-counter+1
        X.sq.S.star[counter]<-stat.star$X.sq.S
        X.sq.S.ij.star[,counter]<-stat.star$X.sq.S.ij
        p.value.ij<-1 - pchisq(q = stat.star$X.sq.S.ij, df = 1)
        p.value.b.min[counter]<-min(p.value.ij)
        p.value.b.prod[counter]<-prod(p.value.ij)
        drop<-0
      }
      # update progress bar
      if (B==B.max) {
        if (print.status) {
          setTxtProgressBar(pb, b)
        }
      }
      if (B!=B.max) {
        if (print.status) {
          left<-B.max-filled
          expect<-(1-drop)*min(left,(B-counter))+drop*max(left,(B.max-b))
          filled<-filled+left/expect
          setTxtProgressBar(pb, filled)
        }
      }
    }
    if (print.status) {
      setTxtProgressBar(pb, B.max)
      close(pb)
    }
    B.use<-min(B,(B.max-discard)) #Only use desired number of resamples
    X.sq.S.star<-X.sq.S.star[1:B.use]
    X.sq.S.ij.star<-X.sq.S.ij.star[,1:B.use]
    p.value.b.min<-p.value.b.min[1:B.use]
    p.value.b.prod<-p.value.b.prod[1:B.use]
    p.value.boot<-mean(X.sq.S.star>=observed$X.sq.S)
    p.combo.min<-list(p = p.value.min, p.star = p.value.b.min, 
                      overall.p = mean(p.value.b.min<=p.value.min))
    p.combo.prod<-list(p = p.value.prod, p.star = p.value.b.prod, 
                       overall.p = mean(p.value.b.prod<=p.value.prod))
    #Histograms
    if (plot.hist) {
      layout.m <- matrix(c(1, 0, 1, 3, 2, 3, 2, 0), nrow = 2, ncol = 4)
      layout(layout.m)
      hist(x = X.sq.S.star, main = "Test using sum statistic", xlab = 
             expression(X[S]^{"2*"}), xlim = c(min(X.sq.S.star, observed$X.sq.S), 
                                               max(X.sq.S.star, observed$X.sq.S)))
      abline(v = observed$X.sq.S, col = "darkgreen", lwd = 5)
      hist(x = p.value.b.prod, main = "Test using product of p-values", xlab = 
             expression(tilde(p)[prod]^{"*"}), xlim = c(min(p.value.b.prod, 
             p.value.prod), max(p.value.b.prod, p.value.prod)))
      abline(v = p.value.prod, col = "darkgreen", lwd = 5)
      hist(x = p.value.b.min, main = "Test using minimum of p-values", xlab = 
             expression(tilde(p)[min]^{"*"}), xlim = c(min(p.value.b.min, 
             p.value.min), max(p.value.b.min, p.value.min)))
      abline(v = p.value.min, col = "darkgreen", lwd = 5) 
    }
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                      summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                      observed$X.sq.S.ij), boot = list(B.use = B.use, 
                      B.discard = discard, p.value.boot = p.value.boot, 
                      p.combo.min.boot = p.combo.min$overall.p, p.combo.prod.boot = 
                      p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star, X.sq.S.ij.star = 
                      X.sq.S.ij.star, p.combo.min.star = p.combo.min$p.star, 
                      p.combo.prod.star = p.combo.prod$p.star))
    output.boot<-list(B.use = B.use, B.discard = discard, 
                      p.value.boot = p.value.boot, p.combo.min.boot = 
                      p.combo.min$overall.p, p.combo.prod.boot = 
                      p.combo.prod$overall.p, X.sq.S.star = X.sq.S.star, 
                      X.sq.S.ij.star = X.sq.S.ij.star, p.combo.min.star = 
                      p.combo.min$p.star, p.combo.prod.star = p.combo.prod$p.star)
  }
  
  #Second-order Rao-Scott adjustment
  if (type == "rs2" | type == "all") {
    #Appendix A calculations
    W.counts<-as.data.frame(table(data[,1:I])) #Get all possible combos of W's
    cols <- c(1:I) #Need to order W's in ascending order starting with W1 
    W.counts.ord<-W.counts[do.call("order", as.data.frame(W.counts[,cols])), ]
    Y.counts<-as.data.frame(table(data[,(I+1):(I+J)])) #All possible Y's
    cols <- c(1:J) #Need to order Y's in ascending order starting with Y1
    Y.counts.ord<-Y.counts[do.call("order", as.data.frame(Y.counts[,cols])), ]
    n.counts<-as.data.frame(table(data)) #Get all possible combos of W's and Y's
    cols <- c(1:(I+J)) #Need to order W's and Y's in ascending order W1, W2, ...
    n.counts.ord<-n.counts[do.call("order", as.data.frame(n.counts[,cols])), ]
    G<-t(data.matrix(W.counts.ord[,1:I])-1) #rx2^r matrix 
    H<-t(data.matrix(Y.counts.ord[,1:J])-1) #cx2^c matrix 
    tau<-n.counts.ord[,(I+J+1)]/n #Vector of multinomial probabilities
    m.row<-G%*%W.counts.ord[,(I+1)] #Vector of W marginal counts
    m.col<-H%*%Y.counts.ord[,(J+1)] #Vector of Y marginal counts
    GH<-kronecker(G,H)  
    m<-GH%*%n.counts.ord[,(I+J+1)] #Marginal table counts (row 1, row 2, ...)
    pi<-m/n   #pi_ij
    pi.row<-m.row/n #pi_i.
    pi.col<-m.col/n #pi_.j  
    j.2r<-matrix(data = 1, nrow = 2^I, ncol = 1) #2^r vector of 1's
    i.2r<-diag(2^I) #(2^r)x(2^r) identity matrix
    j.2c<-matrix(data = 1, nrow = 2^J, ncol = 1) #2^c vector of 1's
    i.2c<-diag(2^J) #(2^c)x(2^c) identity matrix
    G.ij<-G%*%kronecker(i.2r,t(j.2c))
    H.ji<-H%*%kronecker(t(j.2r),i.2c)
    F<-GH - kronecker(pi.row,H.ji) - kronecker(G.ij,pi.col)
    mult.cov<-diag(tau) - tcrossprod(tau)
    sigma<-F%*%tcrossprod(mult.cov,F)
    D<-diag(as.vector(kronecker(pi.row,pi.col)*kronecker(1-pi.row,1-pi.col)))
    Di.sigma<-diag(1/diag(D),dim(D))%*%sigma
    Di.sigma.eigen<-Re(eigen(Di.sigma)$values) #Only use real part of eigenvalues
    sum.Di.sigma.eigen.sq<-sum(Di.sigma.eigen^2)
    X.sq.S.rs2<-I*J*observed$X.sq.S/sum.Di.sigma.eigen.sq 
    df.rs2<-I^2*J^2/sum.Di.sigma.eigen.sq 
    X.sq.S.p.value.rs2<-1-pchisq(q = X.sq.S.rs2, df = df.rs2) 
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                     summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                     observed$X.sq.S.ij), rs2 = list(X.sq.S.rs2 = X.sq.S.rs2, df.rs2 = 
                     df.rs2, p.value.rs2 = X.sq.S.p.value.rs2))
    output.rs2<-list(X.sq.S.rs2 = X.sq.S.rs2, df.rs2 = df.rs2, 
                     p.value.rs2 = X.sq.S.p.value.rs2)
  }
  
  #Bonferroni
  if (type == "bon" | type == "all") {   
    p.value.bon<-min(p.value.min*I*J,1) #p-value should not be greater than 1
    X.sq.S.ij.p.bon<-apply(X = I*J*(1 - pchisq(q = observed$X.sq.S.ij, 
                     df = 1)), MARGIN = c(1,2), FUN = check.min)
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                     summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                     observed$X.sq.S.ij), bon = list(p.value.bon = p.value.bon, 
                     X.sq.S.ij.p.bon = X.sq.S.ij.p.bon))
    output.bon<-list(p.value.bon = p.value.bon, X.sq.S.ij.p.bon = X.sq.S.ij.p.bon) 
  }
  if (type == "all") {
    output<-list(general = list(data = data, I = I, J = J, summary.data = 
                     summary.data, X.sq.S = observed$X.sq.S, X.sq.S.ij = 
                     observed$X.sq.S.ij), boot = output.boot, rs2 = output.rs2, 
                     bon = output.bon)
  }
  class(output)<-"SPMI"
  output
}

#####################################################################################

# MI.test()
# A general function that performs bootstrap testing, the second-order Rao-Scott 
#   approach, and Bonferroni adjustments

MI.test<-function(data, I, J, type = "all", B = 1999, B.max = B, 
                    summary.data = FALSE, add.constant = .5, plot.hist = FALSE, 
                    print.status = TRUE) {
  op<-options()
  on.exit(options(op))
  options(warn=1)
  if((class(data)!="data.frame")&(class(data)!="matrix")) {
    stop("The \"data\" argument requires an object of class \"data.frame\".")
  }
  data<-as.data.frame(data)
  if ((I==0)|(J==0)) {
    stop("\"I\" and \"J\" must be greater than 0.")
  }
  if (!is.numeric(I)|!is.numeric(J)) {
    stop("\"I\" and \"J\" must be numeric.")
  }
  if ((I%%1!=0)|(J%%1!=0)) {
    stop("\"I\" and \"J\" must be integers.")
  }
  I<-as.integer(I)
  J<-as.integer(J)
  if (length(type)>1) {
    warning("The \"type\" argument requires an object of length 1. \n  Only the first element is used.")
    type<-type[1]
  }
  if (!(type%in%c("all","boot","rs2","bon"))) {
    warning("The \"type\" argument can only take on values of \"boot\", \"rs2\", \"bon\", and \"all\". \n  The input value has been changed to the default value of \"all\".")
    type<-"all"
  }
  if (!is.numeric(B)|!is.numeric(B.max)) {
    warning("\"B\" and \"B.max\" must be numeric. \n  The input values have been changed to the default value of 1999.")
    B<-1999
    B.max=B
  }
  if (B.max<B) {
    warning("\"B.max\" must be greater than or equal to \"B\". \n  \"B.max\" has been set equal to \"B\".")
    B.max=B
  }
  if ((B%%1!=0)|(B.max%%1!=0)) {
    warning("\"B\" and \"B.max\" must be integers. \n  The input values have been rounded up to the nearest integer.")
  }
  B<-as.integer(ceiling(B))
  B.max<-as.integer(ceiling(B.max))
  if (!is.numeric(add.constant)) {
    if (add.constant==FALSE) {
      add.constant<-0
    }
    if (add.constant!=FALSE) {
      warning("The \"add.constant\" argument only accepts numeric values. \n  The input value has been changed to the default value of .5.")
      add.constant<-.5
    }
  }
  if (add.constant<0) {
    warning("The \"add.constant\" argument cannot be negative. \n  The input value has been changed to the default value of .5.")
    add.constant<-.5
  }
  if ((class(summary.data)!="logical")&(summary.data!=1)&(summary.data!=0)) {
    stop("The \"summary.data\" argument requires a logical value.")
  }
  if ((class(plot.hist)!="logical")&(plot.hist!=1)&(plot.hist!=0)) {
    warning("The \"plot.hist\" argument requires a logical value. \n  The input value has been changed to the default value of FALSE.")
    plot.hist<-FALSE
  }
  if ((class(print.status)!="logical")&(print.status!=1)&(print.status!=0)) {
    warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    print.status<-TRUE
  }
  summary.data<-as.logical(summary.data)
  plot.hist<-as.logical(plot.hist)
  print.status<-as.logical(print.status)
  #Summary data can only be used with the Bonferroni adjustment
  if (summary.data) { 
    if (type!="bon") {
      warning("Only the Bonferroni adjustment is available when summary.data = TRUE. \n  The \"type\" argument has been changed to \"bon\".")
      type<-"bon"
    }
  }
  nvars<-1+((I>1)&(J>1))
  if (nvars==1) {
    output<-MMI.test(data, I, J, type, B = B, B.max = B.max, 
                     summary.data = summary.data, add.constant = add.constant, plot.hist = plot.hist, 
                     print.status = print.status)
  }
  if (nvars==2) {
    output<-SPMI.test(data, I, J, type, B = B, B.max = B.max, 
                     summary.data = summary.data, add.constant = add.constant, plot.hist = plot.hist, 
                     print.status = print.status)
  }
  output
}

#####################################################################################