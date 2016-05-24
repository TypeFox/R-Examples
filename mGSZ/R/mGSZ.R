# Liisa Holm's Bioinformatics Group, Institute of Biotechnology, University of Helsinki, Finland.
# This software is supplied without any warranty or guarantee.
# Contact: pashupati.mishra@helsinki.fi, petri.toronen@.helsinki.fi

# mGSZ: Modified Gene Set Z-score method for gene set analysis.




mGSZ <- function(x,y,l,f=FALSE,s="T",log=TRUE,g=FALSE,min.sz=5,o=FALSE,pv=0,w1=0.2,w2=0.5,vc=10,p=200){



    expr.data <- x
    gene.sets <- y
    sample.labels <- l
    flip.gene.sets <- f
    select <- s
    is.log <- log
    gene.perm <- g
    min.cl.sz <- min.sz
    other.methods <- o
    pre.var <- pv
    wgt1 <- w1
    wgt2 <- w2
    var.constant <- vc
    start.val=5
    perm.number <- p

   ########################################
   ## Testing that data sizes are correct #
   ########################################
   
   expr.dat.sz <- dim(expr.data)
   gene.sets <- toMatrix(expr.data,gene.sets,flip.gene.sets)
   num.genes <- nrow(gene.sets)
   
   if(!(num.genes == expr.dat.sz[1]) ){
      stop("Number of genes in gene set data and expression data do not match")
   }
   	
    ##Remove rows with no gene set members 
    
	data <- rm.rows.with.noSetMembers(expr.data,gene.sets) 
	expr.data <- data$expr.data
	gene.sets <- data$gene.sets
    
    ##Remove too small gene sets
    
    gene.sets <- rm.small.genesets(gene.sets,min.cl.sz)
    geneset.classes <- colnames(gene.sets)
    num.classes <- dim(gene.sets)[2]
    num.genes <- dim(expr.data)[1]
    geneset.class.sizes <- colSums(gene.sets)
   
    
    
    
    
    ##Calculation of gene scores 


    if(select=="FC"){
        tmp = diffFCscore(expr.data, sample.labels, perm.number)
        tmp.expr.data = tmp$diff_FC
        perm.number = tmp$perm.number
    	}
    else
        tmp = diffScore(expr.data, sample.labels, perm.number)
            if(select=="T"){
                tmp.expr.data = tmp$t_scores}
            if(select=="P"){
                tmp.expr.data=tmp$p_values}
            perm.number = tmp$perm.number
		 
    diff.expr.dat.sz <- dim(tmp.expr.data)
	
	
	## mGSZ score for positive data ####
	
    
	pos.scores <- mGSZ.test.score2(expr.data=tmp.expr.data[,1], gene.sets, wgt1, wgt2, pre.var, var.constant, start.val)
    
	pos.mGSZ.scores <- pos.scores$gene.set.scores$mGSZ.scores

  
    ## Gene set scores for positive data with other methods ##
    
	if(other.methods == TRUE){
		pos.mAllez.scores <- pos.scores$gene.set.scores$mAllez.scores
		pos.mGSA.scores <- pos.scores$gene.set.scores$mGSA.scores
        pos.ss.scores <- SS.test.score(tmp.expr.data[,1], gene.sets)
        pos.sum.scores <- sumTestscore(tmp.expr.data[,1], gene.sets)
        pos.wrs.scores <- WRS.test.score(tmp.expr.data[,1], gene.sets)
        pos.ks.scores <- KS.score(tmp.expr.data[,1], gene.sets)$KS.org
        pos.wks.scores <- KS.score(tmp.expr.data[,1], gene.sets)$KS.mod
        pos.ks.up <- KS.score(tmp.expr.data[,1], gene.sets)$profile2

    }
	 
	
	## mGSZ for sample permutation data ###
	
	
	col.perm.mGSZ <- matrix(0, diff.expr.dat.sz[2]-1, num.classes)
	col.perm.mAllez <- col.perm.mGSZ
	col.perm.mGSA <- col.perm.mGSZ
    col.perm.WRS <- col.perm.mGSZ
    col.perm.SS <- col.perm.mGSZ
    col.perm.SUM <- col.perm.mGSZ
    col.ks <- col.perm.mGSZ
    col.wks <- col.perm.mGSZ
    

    for( k in 1:(diff.expr.dat.sz[2]-1)){
    	tmp <- mGSZ.test.score(tmp.expr.data[,k+1], gene.sets,wgt1, wgt2, pre.var, var.constant, start.val)
    	col.perm.mGSZ[k,]<-tmp$mGSZ.scores

        if(other.methods==TRUE){
    		col.perm.mAllez[k,] <- tmp$mAllez.scores
    		col.perm.mGSA[k,] <- tmp$mGSA.scores
            col.ks[k,] <- KS.score(tmp.expr.data[,k+1], gene.sets)$KS.org
            col.wks[k,] <- KS.score(tmp.expr.data[,k+1], gene.sets)$KS.mod
            col.perm.WRS[k,] <- WRS.test.score(tmp.expr.data[,k+1], gene.sets)
            col.perm.SS[k,] <- SS.test.score(tmp.expr.data[,k+1], gene.sets)
            col.perm.SUM[k,] <- sumTestscore(tmp.expr.data[,k+1], gene.sets)
        }
     }
     
     


### p-value calculation for the gene set scores with sample permutation ###

    mGSZ.p.vals.col.perm <- mGSZ.p.values(pos.mGSZ.scores,col.perm.mGSZ)

    if(other.methods==TRUE){   
        wrs.p.vals.col <- WRS.p.values(pos.wrs.scores,col.perm.WRS)
        ss.p.vals.col <- SS.p.values(pos.ss.scores,col.perm.SS)
        sum.p.vals.col <- mAllez.p.values(pos.sum.scores,col.perm.SUM)
        KS.p.vals.col <- KS.p.values(pos.ks.scores,col.ks)
        wKS.p.vals.col <- KS.p.values(pos.wks.scores,col.wks)
        mGSA.p.vals.col.perm <- mGSA.p.values(pos.mGSA.scores,col.perm.mGSA)
        mAllez.p.vals.col.perm <- mAllez.p.values(pos.mAllez.scores,col.perm.mAllez)
    }



### preparing output table for each of the methods with column permutation ###


    mGSZ.table.col <- data.frame(gene.sets = colnames(gene.sets)[mGSZ.p.vals.col.perm$class.ind], set.size = geneset.class.sizes[mGSZ.p.vals.col.perm$class.ind],gene.set.scores=pos.mGSZ.scores[mGSZ.p.vals.col.perm$class.ind],pvalue=mGSZ.p.vals.col.perm$EV.class,row.names=NULL)



    if(other.methods==TRUE){
        mGSA.table.col <- data.frame(gene.sets = colnames(gene.sets)[mGSA.p.vals.col.perm$class.ind],set.size = geneset.class.sizes[mGSA.p.vals.col.perm$class.ind],gene.set.scores=pos.mGSA.scores[mGSA.p.vals.col.perm$class.ind],pvalue=mGSA.p.vals.col.perm$EV.class,row.names=NULL)
        
   
        mAllez.table.col <- data.frame(gene.sets = colnames(gene.sets)[mAllez.p.vals.col.perm$class.ind],set.size = geneset.class.sizes[mAllez.p.vals.col.perm$class.ind],gene.set.scores=pos.mAllez.scores[mAllez.p.vals.col.perm$class.ind],pvalue=mAllez.p.vals.col.perm$NORM.class,row.names=NULL)
   
        WRS.table.col <- data.frame(gene.sets = colnames(gene.sets)[wrs.p.vals.col$class.ind],set.size = geneset.class.sizes[wrs.p.vals.col$class.ind],gene.set.scores=pos.wrs.scores[wrs.p.vals.col$class.ind],pvalue=wrs.p.vals.col$EMP.class,row.names=NULL)
    
        SUM.table.col <- data.frame(gene.sets = colnames(gene.sets)[sum.p.vals.col$class.ind],set.size = geneset.class.sizes[sum.p.vals.col$class.ind],gene.set.scores=pos.sum.scores[sum.p.vals.col$class.ind],pvalue=sum.p.vals.col$NORM.class,row.names=NULL)
   
        SS.table.col <- data.frame(gene.sets = colnames(gene.sets)[ss.p.vals.col$class.ind],set.size = geneset.class.sizes[ss.p.vals.col$class.ind],gene.set.scores=pos.ss.scores[ss.p.vals.col$class.ind],pvalue=ss.p.vals.col$EMP.class,row.names=NULL)
 
        Old.KS.table.col <- data.frame(gene.sets = colnames(gene.sets)[KS.p.vals.col$class.ind],set.size = geneset.class.sizes[KS.p.vals.col$class.ind],gene.set.scores=pos.ks.scores[KS.p.vals.col$class.ind],pvalue =KS.p.vals.col$EMP.class,row.names=NULL)
    
        New.KS.table.col <- data.frame(gene.sets = colnames(gene.sets)[wKS.p.vals.col$class.ind],set.size = geneset.class.sizes[wKS.p.vals.col$class.ind],gene.set.scores=pos.wks.scores[wKS.p.vals.col$class.ind],pvalue=wKS.p.vals.col$EMP.class,row.names=NULL)
    }
    
    
### Ordering of the tables

    mGSZ.table.col <- mGSZ.table.col[order(mGSZ.table.col$pvalue,decreasing=FALSE),]
    if(other.methods==TRUE){
        mGSA.table.col <- mGSA.table.col[order(mGSA.table.col$pvalue,decreasing=FALSE),]
        mAllez.table.col <- mAllez.table.col[order(mAllez.table.col$pvalue,decreasing=FALSE),]
        WRS.table.col <- WRS.table.col[order(WRS.table.col$pvalue,decreasing=FALSE),]
        SUM.table.col <- SUM.table.col[order(SUM.table.col$pvalue,decreasing=FALSE),]
        SS.table.col <- SS.table.col[order(SS.table.col$pvalue,decreasing=FALSE),]
        Old.KS.table.col <- Old.KS.table.col[order(Old.KS.table.col$pvalue,decreasing=FALSE),]
        New.KS.table.col <- New.KS.table.col[order(New.KS.table.col$pvalue,decreasing=FALSE),]
    }


### calculation of p-value with gene (row) permutation  
    if(gene.perm==TRUE){
        row.perm.mGSZ <- matrix(0, perm.number, num.classes)
        row.perm.mAllez <- row.perm.mGSZ
        row.perm.mGSA <- row.perm.mGSZ
        row.perm.WRS <- row.perm.mGSZ
        row.perm.SS <- row.perm.mGSZ
        row.perm.SUM <- row.perm.mGSZ
        row.ks <- row.perm.mGSZ
        row.wks <- row.perm.mGSZ
        row.permuted.data <- replicate(perm.number,sample(tmp.expr.data[,1], num.genes, replace=FALSE))
        unique.perm <- unique(t(row.permuted.data))
        if(dim(unique.perm)[1]<perm.number){
        #print("Number of unique permutations:" dim(unique.perm)[1])
            row.permuted.data <- t(unique.perm)
        }
        for(i in 1:perm.number){
            res <- mGSZ.test.score4(row.permuted.data[,i],gene.sets,Z_var1=pos.scores$var.attributes$Z_var1, Z_var2=pos.scores$var.attributes$Z_var2, Z_mean1=pos.scores$var.attributes$Z_mean1, Z_mean2=pos.scores$var.attributes$Z_mean2,class_size_index1=pos.scores$var.attributes$class_size_index1, class_size_index2=pos.scores$var.attributes$class_size_index2,start.val)
            row.perm.mGSZ[i,] <- res$mGSZ.scores
            if(other.methods==TRUE){
                row.perm.mAllez[i,] <- res$mAllez.scores
                row.perm.mGSA[i,] <- res$mGSA.scores
                row.ks[k,] <- KS.score(row.permuted.data[,i], gene.sets)$KS.org
                row.wks[k,] <- KS.score(row.permuted.data[,i], gene.sets)$KS.mod
                row.perm.WRS[k,] <- WRS.test.score(row.permuted.data[,i], gene.sets)
                row.perm.SS[k,] <- SS.test.score(row.permuted.data[,i], gene.sets)
                row.perm.SUM[k,] <- sumTestscore(row.permuted.data[,i], gene.sets)
            }

		}
	

    ### p-value calculation for the gene set scores with gene (row) permutation ###

        mGSZ.p.vals.row.perm <- mGSZ.p.values(pos.mGSZ.scores,row.perm.mGSZ)
        if(other.methods==TRUE){
            wrs.p.vals.row <- WRS.p.values(pos.wrs.scores,row.perm.WRS)
            ss.p.vals.row <- SS.p.values(pos.ss.scores,row.perm.SS)
            sum.p.vals.row <- mAllez.p.values(pos.sum.scores,row.perm.SUM)
            KS.p.vals.row <- KS.p.values(pos.ks.scores,row.ks)
            wKS.p.vals.row <- KS.p.values(pos.wks.scores,row.wks)
            mGSA.p.vals.row.perm <- mGSA.p.values(pos.mGSA.scores,row.perm.mGSA)
            mAllez.p.vals.row.perm <- mAllez.p.values(pos.mAllez.scores,row.perm.mAllez)
        }



        ##### preparing output table for each of the methods with row permutation ###
        
        mGSZ.table.row <- data.frame(gene.sets = colnames(gene.sets)[mGSZ.p.vals.row.perm$class.ind],set.size = geneset.class.sizes[mGSZ.p.vals.row.perm$class.ind],gene.set.scores=pos.mGSZ.scores[mGSZ.p.vals.row.perm$class.ind],pvalue=mGSZ.p.vals.row.perm$EV.class,row.names=NULL)
        
    
        if(other.methods==TRUE){
            mGSA.table.row <- data.frame(gene.sets = colnames(gene.sets)[mGSA.p.vals.row.perm$class.ind],set.size = geneset.class.sizes[mGSA.p.vals.row.perm$class.ind],gene.set.scores=pos.mGSA.scores[mGSA.p.vals.row.perm$class.ind],pvalue=mGSA.p.vals.row.perm$EV.class,row.names=NULL)
      
            mAllez.table.row <- data.frame(gene.sets = colnames(gene.sets)[mAllez.p.vals.row.perm$class.ind],set.size = geneset.class.sizes[mAllez.p.vals.row.perm$class.ind],gene.set.scores=pos.mAllez.scores[mAllez.p.vals.row.perm$class.ind],pvalue=mAllez.p.vals.row.perm$NORM.class,row.names=NULL)
        
            WRS.table.row <- data.frame(gene.sets = colnames(gene.sets)[wrs.p.vals.row$class.ind],set.size = geneset.class.sizes[wrs.p.vals.row$class.ind],gene.set.scores=pos.wrs.scores[wrs.p.vals.row$class.ind],pvalue=wrs.p.vals.row$EMP.class,row.names=NULL)
        
        
            SUM.table.row <- data.frame(gene.sets = colnames(gene.sets)[sum.p.vals.row$class.ind],set.size = geneset.class.sizes[sum.p.vals.row$class.ind],gene.set.scores=pos.sum.scores[sum.p.vals.row$class.ind],pvalue=sum.p.vals.row$NORM.class,row.names=NULL)
       
            SS.table.row <- data.frame(gene.sets = colnames(gene.sets)[ss.p.vals.row$class.ind],set.size = geneset.class.sizes[ss.p.vals.row$class.ind],gene.set.scores=pos.ss.scores[ss.p.vals.row$class.ind],pvalue=ss.p.vals.row$EMP.class,row.names=NULL)
       
            Old.KS.table.row <- data.frame(gene.sets = colnames(gene.sets)[KS.p.vals.row$class.ind],set.size = geneset.class.sizes[KS.p.vals.row$class.ind],gene.set.scores=pos.ks.scores[KS.p.vals.row$class.ind],pvalue=KS.p.vals.row$EMP.class,row.names=NULL)
    
            New.KS.table.row <- data.frame(gene.sets = colnames(gene.sets)[wKS.p.vals.row$class.ind],set.size = geneset.class.sizes[wKS.p.vals.row$class.ind],gene.set.scores=pos.wks.scores[wKS.p.vals.row$class.ind],pvalue=wKS.p.vals.row$EMP.class,row.names=NULL)
        }
    
    
    ## Ordering of the tables for row permutation

        mGSZ.table.row <- mGSZ.table.row[order(mGSZ.table.row$pvalue,decreasing=FALSE),]
        if(other.methods==TRUE){
            mGSA.table.row <- mGSA.table.row[order(mGSA.table.row$pvalue,decreasing=FALSE),]
            mAllez.table.row <- mAllez.table.row[order(mAllez.table.row$pvalue,decreasing=FALSE),]
            WRS.table.row <- WRS.table.row[order(WRS.table.row$pvalue,decreasing=FALSE),]
            SUM.table.row <- SUM.table.row[order(SUM.table.row$pvalue,decreasing=FALSE),]
            SS.table.row <- SS.table.row[order(SS.table.row$pvalue,decreasing=FALSE),]
            Old.KS.table.row <- Old.KS.table.row[order(Old.KS.table.row$pvalue,decreasing=FALSE),]
            New.KS.table.row <- New.KS.table.row[order(New.KS.table.row$pvalue,decreasing=FALSE),]
        }
        
        if(other.methods==TRUE){
            out <- list(sample.perm =list(mGSZ = mGSZ.table.col,mGSA = mGSA.table.col,mAllez = mAllez.table.col,WRS = WRS.table.col,KS = Old.KS.table.col,wKS = New.KS.table.col,SS = SS.table.col,SUM = SUM.table.col),gene.perm=list(mGSZ = mGSZ.table.row,mGSA = mGSA.table.row,mAllez = mAllez.table.row,WRS = WRS.table.row,KS = Old.KS.table.row,wKS = New.KS.table.row, SS = SS.table.row,SUM = SUM.table.row), sample.labels = sample.labels, perm.number = perm.number, expr.data = expr.data, gene.sets = gene.sets, flip.gene.sets = flip.gene.sets, min.cl.sz = min.cl.sz, other.methods = other.methods, pre.var = pre.var, wgt2 = wgt2, wgt1 = wgt1, var.constant = var.constant, start.val = start.val, select = select, is.log = is.log, gene.perm.log = gene.perm)}
        else{
            out <- list(mGSZ.sample.perm = mGSZ.table.col, mGSZ.gene.perm = mGSZ.table.row, sample.labels = sample.labels, perm.number = perm.number, expr.data = expr.data, gene.sets = gene.sets, flip.gene.sets = flip.gene.sets, min.cl.sz = min.cl.sz, other.methods = other.methods, pre.var = pre.var, wgt2 = wgt2, wgt1 = wgt1, var.constant = var.constant, start.val = start.val, select = select, is.log = is.log, gene.perm.log = gene.perm)
            }
        return(out)


    }
    else{
        if(other.methods==TRUE){
            out <- list(mGSZ = mGSZ.table.col, mGSA = mGSA.table.col, mAllez = mAllez.table.col,WRS = WRS.table.col,KS = Old.KS.table.col,wKS = New.KS.table.col,SUM = SUM.table.col,SS = SS.table.col, sample.labels = sample.labels, perm.number = perm.number, expr.data = expr.data, gene.sets = gene.sets, flip.gene.sets = flip.gene.sets, min.cl.sz = min.cl.sz, other.methods = other.methods, pre.var = pre.var, wgt2 = wgt2, wgt1 = wgt1, var.constant = var.constant, start.val = start.val, select = select, is.log = is.log, gene.perm.log = gene.perm)
                        }
            else{
            out <- list(mGSZ = mGSZ.table.col, sample.labels = sample.labels, perm.number = perm.number, expr.data = expr.data, gene.sets = gene.sets, flip.gene.sets = flip.gene.sets, min.cl.sz = min.cl.sz, other.methods = other.methods, pre.var = pre.var, wgt2 = wgt2, wgt1 = wgt1, var.constant = var.constant, start.val = start.val, select = select, is.log = is.log, gene.perm.log = gene.perm)
            }
        return(out)
    }

}




################

rm.rows.with.noSetMembers <-
function(expr.data,gene.sets){
	tmp_sum <- apply(gene.sets, 1, sum)
    ind_wth_sets <- which(tmp_sum > 0)
    expr.data <- expr.data[ind_wth_sets,]
    gene.sets   <- gene.sets[ind_wth_sets,]
    out <- list(expr.data=expr.data, gene.sets=gene.sets)}
    
###############

rm.small.genesets <-
function(gene.sets,min.cl.sz){
    tmp_sum <- apply(gene.sets, 2, sum)
    ind_wth_sets <- which(tmp_sum > min.cl.sz)
    gene.sets   <- gene.sets[,ind_wth_sets]
	out=gene.sets}

##################

sumVarMean_calc <-
function(expr_data, gene.sets, pre.var){
	dim_sets <- dim(gene.sets)
	length_expr <- length(expr_data)
	set_sz <- colSums(gene.sets)

	unique_class_sz <- unique(set_sz)
	num_genes <- length(expr_data)
	
	divider <- c(1:num_genes)
    mean_table <- cumsum(expr_data)/divider
    mean_table_sq <- mean_table^2
    var_table <- cumsum(expr_data^2)/divider - (mean_table)^2 + pre.var
    class_sz_index <- rep(0,dim_sets[2])
    for (i in 1:dim_sets[2]){
    	class_sz_index[i] <- which(set_sz[i]==unique_class_sz)
    }
    var_table <- var_table + pre.var
    hyge_stat <- count_hyge_var_mean(dim_sets[1],unique_class_sz)
    
    z_var <- matrix(numeric(0),dim_sets[1],length(unique_class_sz)) 
    
    z_mean <- z_var
    max_value <- 1:dim_sets[1]
    for (j in 1:length(unique_class_sz)){
    	prob_sum <- count.prob.sum(dim_sets[1],hyge_stat$mean[,j],hyge_stat$var[,j])
    	z_var[,j] <- 4*(var_table*prob_sum + mean_table_sq*hyge_stat$var[,j])
    	z_mean[,j] <- mean_table*(2*hyge_stat$mean[,j]-max_value)}
        
        out = list(Z_var = z_var, Z_mean = z_mean, class_size_index = class_sz_index,var_table=var_table,mean_table_sq=mean_table_sq,set_sz=set_sz)
        return(out)
  
    
}

#########

count.prob.sum <-
function(M, hyge.mean, hyge.var){
    tulos = matrix(rep(0,length(hyge.mean)),byrow=FALSE)
    N_tab = matrix(c(2:M),byrow=FALSE)
    tulos[1] = hyge.mean[1]
    tulos[2:M] = N_tab/(N_tab-1)*hyge.mean[2:M]-(hyge.mean[2:M]^2+hyge.var[2:M])/(N_tab-1)
    tulos = tulos
}

############

count_hyge_var_mean <-
function(M,K){
    N = matrix(rep((1:M),length(K)),ncol=length(K))
    K = matrix(rep(K,M),byrow=TRUE,ncol=length(K))
    out1 = N*K/M
    out2 = N*(K*(M-K)/M^2)*((M-N)/(M-1))
    results = list(mean= out1, var= out2)
}

##############

calc_z_var <-
function(num.genes,unique_class_sz_ln,pre_z_var,wgt2,var.constant){
	ones_matrix = matrix(1,num.genes,unique_class_sz_ln)
	ones_tmatrix = t(ones_matrix)
	median_matrix = apply(pre_z_var,2,median)
	pre_median_part = ones_tmatrix*median_matrix*wgt2
	median_part = t(pre_median_part)
	z_var = (pre_z_var + median_part + var.constant)^0.5
}

#############

geneSetsList <-
function(data){
  data <- GSA.read.gmt(data)
  geneSets <- data$genesets
  names(geneSets) <- data$geneset.names
  return(geneSets)
}






