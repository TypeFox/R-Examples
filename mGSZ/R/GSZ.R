mGSZ.test.score <-
function(expr.data, gene.sets, wgt1 = 0.2, wgt2 = 0.5, pre.var = 0, var.constant = 10, start.val = 5, flip.gene.sets = FALSE, table = FALSE,...){

    num.genes <- length(expr.data)

  # Ordering of gene expression data and gene sets data

  gene.sets <- toMatrix(expr.data,gene.sets,flip.gene.sets)
  ord_out <- order(expr.data, decreasing= TRUE)
  expr.data <- expr.data[ord_out]
  set.dim <- dim(gene.sets)
  cols <- set.dim[2]
  gene.sets <- gene.sets[ord_out,]
 
    
  expr.data.ud <- expr.data[num.genes:1] # expression values turned up-side down
                                         # This does the analysis of the lower end


	num.genes=length(expr.data)
	num.classes=dim(gene.sets)[2]
	set_sz <- apply(gene.sets,2,sum)
	unique_class_sz_ln <- length(unique(set_sz))
	
	
	### Defining variables for different output

	mAllez <- rep(0,num.classes)
	mgsa <- matrix(0,2,num.classes)

	pre_z_var.1 <- sumVarMean_calc(expr.data, gene.sets, pre.var)
	pre_z_var.2 <- sumVarMean_calc(expr.data.ud, gene.sets, pre.var)
	
	Z_var1 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.1$Z_var,wgt2,var.constant)
	Z_var2 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.2$Z_var,wgt2,var.constant)

	out2 = matrix(0,2,num.classes)
	mgsa1 <- numeric(num.classes)
	mgsa2 <- mgsa1
	mgsz1 <- mgsa1
	mgsz2 <- mgsa2
	for (k in 1:num.classes){
		po1 <- which(gene.sets[,k]==1)
		po0 <- which(gene.sets[,k]==0)
		if(length(po1) > 0 & length(po0) > 0){
			tmp1 = expr.data
			tmp1[po0] = 0
			tmp0 = expr.data
			tmp0[po1] = 0
			result1 = cumsum(tmp1)-cumsum(tmp0)-pre_z_var.1$Z_mean[,pre_z_var.1$class_size_index[k]]
			result2 = cumsum(tmp1[num.genes:1])-cumsum(tmp0[num.genes:1])-pre_z_var.2$Z_mean[,pre_z_var.2$class_size_index[k]]
			
			result1[1:start.val] <- 0
			result2[1:start.val] <- 0
			A = result1/Z_var1[,pre_z_var.1$class_size_index[k]]
			B = result2/Z_var2[,pre_z_var.2$class_size_index[k]]
			mAllez[k] <- A[num.genes]
			mgsa1[k] <- A[round(num.genes/2)]
			mgsa2[k] <- B[round(num.genes/2)]
			mgsz1[k] <- max(abs(A))
			mgsz2[k] <- max(abs(B))
			
	    }
    }
    mgsa[1,]<-mgsa1
    mgsa[2,]<-mgsa2
    out2[1,]<-mgsz1
    out2[2,]<-mgsz2
    if(table){
        out <- cbind(Gene.sets = colnames(gene.sets), mGSZ.scores = apply(out2,2,max), mAllez.scores = mAllez, mGSA.scores = apply(abs(mgsa),2,max))
    }
    else{
    out = list(mGSZ.scores = apply(out2,2,max), mAllez.scores = mAllez, mGSA.scores= apply(abs(mgsa),2,max))
    return(out)}
    
 }
 
 ####################
 
 mGSZ.test.score.stb2 <-
function(gene.sets.stb, expr.data, gene.sets, wgt1, wgt2, pre.var, var.constant,start.val){
	

  num.genes <- length(expr.data)

  # Ordering of gene expression data and gene sets

  ord_out <- order(expr.data, decreasing= TRUE)
  expr.data <- expr.data[ord_out]
  set.dim <- dim(gene.sets)
  cols <- set.dim[2]
  gene.sets <- gene.sets[ord_out,] 
                           
 
  expr.data.ud <- expr.data[num.genes:1] 
                                         


	num.genes=length(expr.data)
	num.classes=dim(gene.sets)[2]
	set_sz <- apply(gene.sets,2,sum)
	unique_class_sz_ln <- length(unique(set_sz))
	
	
	### Defining variables for different output
	
    mGSZ.up <- list()
    #GSZ.down <- list()

    pre_z_var.1 <- sumVarMean_calc(expr.data, gene.sets, pre.var)
	pre_z_var.2 <- sumVarMean_calc(expr.data.ud, gene.sets, pre.var)
	
	Z_var1 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.1$Z_var,wgt2,var.constant)
	Z_var2 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.2$Z_var,wgt2,var.constant)

	out2 = matrix(0,2,num.classes)
	mgsz1 <- numeric(num.classes)
	mgsz2 <- mgsz1
    
	
	for (k in 1:num.classes){
		po1 <- which(gene.sets[,k]==1)
		po0 <- which(gene.sets[,k]==0)
		if(length(po1) > 0 & length(po0) > 0){
			tmp1 = expr.data
			tmp1[po0] = 0
			tmp0 = expr.data
			tmp0[po1] = 0
			result1 = cumsum(tmp1)-cumsum(tmp0)-pre_z_var.1$Z_mean[,pre_z_var.1$class_size_index[k]]
			result2 = cumsum(tmp1[num.genes:1])-cumsum(tmp0[num.genes:1])-pre_z_var.2$Z_mean[,pre_z_var.2$class_size_index[k]]
			
			result1[1:start.val] <- 0
			result2[1:start.val] <- 0
			A = result1/Z_var1[,pre_z_var.1$class_size_index[k]]
			B = result2/Z_var2[,pre_z_var.2$class_size_index[k]]

            if(colnames(gene.sets)[k] %in% gene.sets.stb){
            
            mGSZ.up[[length(mGSZ.up)+1]] <- (A)
            #GSZ.down[[length(mGSZ.down)+1]] <- (B)
            }
    
			mgsz1[k] <- max(abs(A))
			mgsz2[k] <- max(abs(B))
			
	    }
    }
    
    out2[1,]<-mgsz1
    out2[2,]<-mgsz2
    result1 = list(GSZ.result = apply(out2,2,max))
    result2 = list(Z_var1=Z_var1,Z_var2=Z_var2,Z_mean1=pre_z_var.1$Z_mean,Z_mean2=pre_z_var.2$Z_mean,class_size_index1=pre_z_var.1$class_size_index,class_size_index2=pre_z_var.2$class_size_index)
    
    out = list(gene.set.scores=result1, var.attributes=result2,mGSZ.up=mGSZ.up)

    return(out)
 }

################

mGSZ.test.score.stb3 <-
function(gene.sets.stb, expr.data, gene.sets, wgt1, wgt2, pre.var, var.constant,start.val,...){
	

  num.genes <- length(expr.data)

  # Ordering of gene expression data and gene set data

 
    
      ord_out <- order(expr.data, decreasing= TRUE)
      expr.data <- expr.data[ord_out]
      set.dim <- dim(gene.sets)
      cols <- set.dim[2]
      gene.sets <- gene.sets[ord_out,]
                                          
  

  expr.data.ud <- expr.data[num.genes:1] 
                                         


	num.genes=length(expr.data)
	num.classes=dim(gene.sets)[2]
	set_sz <- apply(gene.sets,2,sum)
	unique_class_sz_ln <- length(unique(set_sz))
	
	
	### Defining variables for different output

    mGSZ.up <- list()
    #GSZ.down <- list()
	
	pre_z_var.1 <- sumVarMean_calc(expr.data, gene.sets, pre.var)
	pre_z_var.2 <- sumVarMean_calc(expr.data.ud, gene.sets, pre.var)
	
	Z_var1 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.1$Z_var,wgt2,var.constant)
	Z_var2 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.2$Z_var,wgt2,var.constant)

	out2 = matrix(0,2,num.classes)
	mgsz1 <- numeric(num.classes)
	mgsz2 <- mgsz1
   

	for (k in 1:num.classes){
		po1 <- which(gene.sets[,k]==1)
		po0 <- which(gene.sets[,k]==0)
		if(length(po1) > 0 & length(po0) > 0){
			tmp1 = expr.data
			tmp1[po0] = 0
			tmp0 = expr.data
			tmp0[po1] = 0
			result1 = cumsum(tmp1)-cumsum(tmp0)-pre_z_var.1$Z_mean[,pre_z_var.1$class_size_index[k]]
			result2 = cumsum(tmp1[num.genes:1])-cumsum(tmp0[num.genes:1])-pre_z_var.2$Z_mean[,pre_z_var.2$class_size_index[k]]
			
			result1[1:start.val] <- 0
			result2[1:start.val] <- 0
			A = result1/Z_var1[,pre_z_var.1$class_size_index[k]]
			B = result2/Z_var2[,pre_z_var.2$class_size_index[k]]
            
            if(colnames(gene.sets)[k] %in% gene.sets.stb){
                
                mGSZ.up[[length(mGSZ.up)+1]] <- (A)
                #mGSZ.down[[length(mGSZ.down)+1]] <- (B)
                }
			mgsz1[k] <- max(abs(A))
            mgsz2[k] <- max(abs(B))
			
	    }
    }
    out2[1,]<-mgsz1
    out2[2,]<-mgsz2
    out = list(mGSZ.result = apply(out2,2,max),mGSZ.up=mGSZ.up)

    return(out)
    
 }

#################

mGSZ.test.score.stb4 <-
function(gene.sets.stb, expr.data,gene.sets,Z_var1,Z_var2,Z_mean1,Z_mean2,class_size_index1,class_size_index2,start.val){


  num.genes <- length(expr.data)

  # Ordering of gene expression data and gene sets 

    
      ord_out <- order(expr.data, decreasing= TRUE)
      expr.data <- expr.data[ord_out]
      set.dim <- dim(gene.sets)
      cols <- set.dim[2]
      gene.sets <- gene.sets[ord_out,] 
                                          
  
  expr.data.ud <- expr.data[num.genes:1] 


	num.genes=length(expr.data)
	num.classes=dim(gene.sets)[2]
	set_sz <- apply(gene.sets,2,sum)
	unique_class_sz_ln <- length(unique(set_sz))
	
	
	### Defining variables for different output

    mGSZ.up <- list()
    #mGSZ.down <- list()
   	out2 = matrix(0,2,num.classes)
	mgsz1 <- numeric(num.classes)
	mgsz2 <- mgsz1

	for (k in 1:num.classes){
		po1 <- which(gene.sets[,k]==1)
		po0 <- which(gene.sets[,k]==0)
		if(length(po1) > 0 & length(po0) > 0){
			tmp1 = expr.data
			tmp1[po0] = 0
			tmp0 = expr.data
			tmp0[po1] = 0
			result1 = cumsum(tmp1)-cumsum(tmp0)-Z_mean1[,class_size_index1[k]]
			result2 = cumsum(tmp1[num.genes:1])-cumsum(tmp0[num.genes:1])-Z_mean2[,class_size_index2[k]]
			
			result1[1:start.val] <- 0
			result2[1:start.val] <- 0
			A = result1/Z_var1[,class_size_index1[k]]
			B = result2/Z_var2[,class_size_index2[k]]
            if(colnames(gene.sets)[k] %in% gene.sets.stb){
            
                mGSZ.up[[length(mGSZ.up)+1]] <- (A)
                #mGSZ.down[[length(mGSZ.down)]] <- (B)
                }
			mgsz1[k] <- max(abs(A))
			mgsz2[k] <- max(abs(B))
			
	    }
    }
    out2[1,]<-mgsz1
    out2[2,]<-mgsz2
    
    out = list(mGSZ.result = apply(out2,2,max),mGSZ.up=mGSZ.up)
    return(out)
}

###############

mGSZ.test.score2 <-
function(expr.data, gene.sets, wgt1, wgt2, pre.var, var.constant,start.val, flip.gene.sets = FALSE,...){
	
  num.genes <- length(expr.data)

  # Ordering of gene expression data and gene sets data

  gene.sets <- toMatrix(expr.data,gene.sets,flip.gene.sets)
  ord_out <- order(expr.data, decreasing= TRUE)
  expr.data <- expr.data[ord_out]
  set.dim <- dim(gene.sets)
  cols <- set.dim[2]
  gene.sets <- gene.sets[ord_out,] 
                           

 
  expr.data.ud <- expr.data[num.genes:1] # expression values turned up-side down
                                         # This does the analysis of the lower end


	num.genes=length(expr.data)
	num.classes=dim(gene.sets)[2]
	set_sz <- apply(gene.sets,2,sum)
	unique_class_sz_ln <- length(unique(set_sz))
	
	mAllez <- rep(0,num.classes)
    
	mgsa <- matrix(0,2,num.classes)

	pre_z_var.1 <- sumVarMean_calc(expr.data, gene.sets, pre.var)
	pre_z_var.2 <- sumVarMean_calc(expr.data.ud, gene.sets, pre.var)
	
	Z_var1 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.1$Z_var,wgt2,var.constant)
	Z_var2 = calc_z_var(num.genes,unique_class_sz_ln, pre_z_var.2$Z_var,wgt2,var.constant)

	out2 = matrix(0,2,num.classes)
	mgsa1 <- numeric(num.classes)
	mgsa2 <- mgsa1
	mgsz1 <- mgsa1
	mgsz2 <- mgsa2
	for (k in 1:num.classes){
		po1 <- which(gene.sets[,k]==1)
		po0 <- which(gene.sets[,k]==0)
		if(length(po1) > 0 & length(po0) > 0){
			tmp1 = expr.data
			tmp1[po0] = 0
			tmp0 = expr.data
			tmp0[po1] = 0
			result1 = cumsum(tmp1)-cumsum(tmp0)-pre_z_var.1$Z_mean[,pre_z_var.1$class_size_index[k]]
			result2 = cumsum(tmp1[num.genes:1])-cumsum(tmp0[num.genes:1])-pre_z_var.2$Z_mean[,pre_z_var.2$class_size_index[k]]
			
			result1[1:start.val] <- 0
			result2[1:start.val] <- 0
			A = result1/Z_var1[,pre_z_var.1$class_size_index[k]]
			B = result2/Z_var2[,pre_z_var.2$class_size_index[k]]
           
            mAllez[k] <- A[num.genes]
			mgsa1[k] <- A[round(num.genes/2)]
			mgsa2[k] <- B[round(num.genes/2)]
			mgsz1[k] <- max(abs(A))
			mgsz2[k] <- max(abs(B))
			
	    }
    }
    mgsa[1,]<-mgsa1
    mgsa[2,]<-mgsa2
    out2[1,]<-mgsz1
    out2[2,]<-mgsz2
    result1 = list(mGSZ.scores = apply(out2,2,max), mAllez.scores = mAllez, mGSA.scores= apply(abs(mgsa),2,max))
    result2 = list(Z_var1=Z_var1,Z_var2=Z_var2,Z_mean1=pre_z_var.1$Z_mean,Z_mean2=pre_z_var.2$Z_mean,class_size_index1=pre_z_var.1$class_size_index,class_size_index2=pre_z_var.2$class_size_index)
    
    out = list(gene.set.scores=result1, var.attributes=result2)
    return(out)
 }
 

################

mGSZ.test.score4 <-
function(expr.data,gene.sets,Z_var1,Z_var2,Z_mean1,Z_mean2,class_size_index1,class_size_index2,start.val,flip.gene.sets=FALSE){

    num.genes <- length(expr.data)

  # Ordering of gene expression data and gene sets data

 gene.sets <- toMatrix(expr.data,gene.sets,flip.gene.sets)
 ord_out <- order(expr.data, decreasing= TRUE)
 expr.data <- expr.data[ord_out]
 set.dim <- dim(gene.sets)
 cols <- set.dim[2]
 gene.sets <- gene.sets[ord_out,]

  expr.data.ud <- expr.data[num.genes:1] # expression values turned up-side down
                                         # This does the analysis of the lower end


	num.genes=length(expr.data)
	num.classes=dim(gene.sets)[2]
	set_sz <- apply(gene.sets,2,sum)
	unique_class_sz_ln <- length(unique(set_sz))
	
	
	### Defining variables for different output
	
	mAllez <- rep(0,num.classes)
	mgsa <- matrix(0,2,num.classes)

	out2 = matrix(0,2,num.classes)
	mgsa1 <- numeric(num.classes)
	mgsa2 <- mgsa1
	mgsz1 <- mgsa1
	mgsz2 <- mgsa2

	for (k in 1:num.classes){
		po1 <- which(gene.sets[,k]==1)
		po0 <- which(gene.sets[,k]==0)
		if(length(po1) > 0 & length(po0) > 0){
			tmp1 = expr.data
			tmp1[po0] = 0
			tmp0 = expr.data
			tmp0[po1] = 0
			result1 = cumsum(tmp1)-cumsum(tmp0)-Z_mean1[,class_size_index1[k]]
			
			result2 = cumsum(tmp1[num.genes:1])-cumsum(tmp0[num.genes:1])-Z_mean2[,class_size_index2[k]]
			
			result1[1:start.val] <- 0
			result2[1:start.val] <- 0
			A = result1/Z_var1[,class_size_index1[k]]
			B = result2/Z_var2[,class_size_index2[k]]
			mAllez[k] <- A[num.genes]
			mgsa1[k] <- A[round(num.genes/2)]
			mgsa2[k] <- B[round(num.genes/2)]
			mgsz1[k] <- max(abs(A))
			mgsz2[k] <- max(abs(B))
			
	    }
    }
    mgsa[1,]<-mgsa1
    mgsa[2,]<-mgsa2
    out2[1,]<-mgsz1
    out2[2,]<-mgsz2
    out = list(mGSZ.scores = apply(out2,2,max), mAllez.scores = mAllez, mGSA.scores= apply(abs(mgsa),2,max))

    return(out)
}

