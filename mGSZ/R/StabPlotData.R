StabPlotData <-
function(mGSZobj,rank.vector,sample.perm.data=FALSE){

# sample.perm.data: TRUE if sample perm data is to be plotted instead of gene data
   ########################################
   ## Testing that data sizes are correct #
   ########################################
   
   gene.sets <- mGSZobj$gene.sets
   expr.data <- mGSZobj$expr.data
   expr.dat.sz <- dim(expr.data)
   perm.number <- mGSZobj$perm.number
   sample.labels <- mGSZobj$sample.labels
   flip.gene.sets <- mGSZobj$flip.gene.sets
   select <- mGSZobj$select
   is.log <- mGSZobj$is.log
   gene.perm <- mGSZobj$gene.perm.log
   min.cl.sz <- mGSZobj$min.cl.sz
   other.methods <- mGSZobj$other.methods
   pre.var <- mGSZobj$pre.var
   wgt1 <- mGSZobj$wgt1
   wgt2 <- mGSZobj$wgt2
   var.constant <- mGSZobj$var.constant
   start.val <- mGSZobj$start.val
   
    
   geneset.classes <- colnames(gene.sets)
   num.classes <- dim(gene.sets)[2]
   num.genes <- dim(expr.data)[1]
   geneset.class.sizes <- colSums(gene.sets)


    #require(limma)
    #require(Biobase)
    if(select=="FC"){
        tmp = diffFCscore(expr.data, sample.labels, perm.number)
        tmp.expr.data = tmp$diff_FC
    
    }
        
      else
            tmp = diffScore(expr.data, sample.labels, perm.number)
                if(select=="T"){
                    tmp.expr.data = tmp$t_scores}
                if(select=="P"){
                    tmp.expr.data=tmp$p_values}
		 
        
	diff.expr.dat.sz <- dim(tmp.expr.data)
    
     if(gene.perm==TRUE & other.methods==TRUE){
        gene.sets.stb <- as.vector(mGSZobj$gene.perm$mGSZ[,1][rank.vector])
        
    }
    
    if(gene.perm==TRUE & other.methods==FALSE){
        gene.sets.stb <- as.vector(mGSZobj$mGSZ.gene.perm[,1][rank.vector])
        
    }

    if(gene.perm==TRUE & sample.perm.data==TRUE & other.methods==TRUE){
        gene.sets.stb <- as.vector(mGSZobj$sample.perm$mGSZ[,1][rank.vector])
    }
    
    if(gene.perm==TRUE & sample.perm.data==TRUE & other.methods==FALSE){
        gene.sets.stb <- as.vector(mGSZobj$mGSZ.sample.perm[,1][rank.vector])
    }
    
    if(gene.perm==FALSE & other.methods==TRUE){
        gene.sets.stb <- as.vector(mGSZobj$mGSZ[,1][rank.vector])
    }
    
    if(gene.perm==FALSE & other.methods==FALSE){
        gene.sets.stb <- as.vector(mGSZobj$mGSZ[,1][rank.vector])
    }

	
	## mGSZ for positive data ####
	#############################
	pos.scores <- mGSZ.test.score.stb2(gene.sets.stb=gene.sets.stb,expr.data=tmp.expr.data[,1], gene.sets,wgt1=wgt1,wgt2=wgt2,pre.var=pre.var,var.constant=var.constant,start.val=start.val)
    mGSZ.up <- pos.scores$mGSZ.up

    ## mGSZ for col randomized data ###
	##################################
	
	col.perm.mGSZ <- matrix(0, diff.expr.dat.sz[2]-1, num.classes)
    perm.number <- diff.expr.dat.sz[2]-1
    mGSZ.cperm.up <- vector("list",length(rank.vector))

	
    for( k in 1:(diff.expr.dat.sz[2]-1)){
    	tmp <- mGSZ.test.score.stb3(gene.sets.stb=gene.sets.stb,expr.data=tmp.expr.data[,k+1], gene.sets=gene.sets,wgt1=wgt1,wgt2=wgt2,pre.var=pre.var,var.constant=var.constant,start.val=start.val)
        mGSZ.cperm.up[[k]] <- tmp$mGSZ.up

     }
     
     out <- list(positive.prof.up = mGSZ.up, perm.data = mGSZ.cperm.up, expr.data = expr.data, gene.sets = gene.sets, perm.number = perm.number,rank.vector=rank.vector,gene.sets.stb=gene.sets.stb)
     
     

    if(gene.perm==TRUE){
        mGSZ.rperm.up <- vector("list",length(rank.vector))
        row.permuted.data <- replicate(perm.number,sample(tmp.expr.data[,1], num.genes, replace=FALSE))
        unique.perm <- unique(t(row.permuted.data))
        row.permuted.data <- t(unique.perm)
        
        
        for(i in 1:perm.number){
            res <- mGSZ.test.score.stb4(gene.sets.stb=gene.sets.stb, expr.data=row.permuted.data[,i],gene.sets=gene.sets,Z_var1=pos.scores$var.attributes$Z_var1, Z_var2=pos.scores$var.attributes$Z_var2, Z_mean1=pos.scores$var.attributes$Z_mean1, Z_mean2=pos.scores$var.attributes$Z_mean2,class_size_index1=pos.scores$var.attributes$class_size_index1, class_size_index2=pos.scores$var.attributes$class_size_index2,start.val)
            mGSZ.rperm.up[[i]] <- res$mGSZ.up
            }
            
            out <- list(positive.prof.up = mGSZ.up, perm.data = mGSZ.rperm.up, expr.data = expr.data, gene.sets = gene.sets, perm.number = perm.number,rank.vector=rank.vector,gene.sets.stb=gene.sets.stb)
            if(gene.perm==TRUE & sample.perm.data==TRUE){
                out <- list(positive.prof.up = mGSZ.up, perm.data = mGSZ.cperm.up, expr.data = expr.data, gene.sets = gene.sets, perm.number = perm.number,rank.vector=rank.vector,gene.sets.stb=gene.sets.stb)
            }
	}
out
}
