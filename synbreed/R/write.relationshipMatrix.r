write.relationshipMatrix <- function(x,file=NULL,sorting=c("WOMBAT","ASReml"),type=c("ginv","inv","none"),digits=10){
        
        type <- match.arg(type)
        sorting <- match.arg(sorting)

        if(sorting=="WOMBAT" & type!="ginv") stop("'type' must be 'ginv' for WOMBAT")
        
        # pass (inverse) relationship matrix
        if(type=="ginv") rMinv <- ginv(x)
        if(type=="inv")  rMinv <- solve(x)
        if(type=="none") rMinv <- x
        
        rMinv <- round(rMinv,digits)
        
        # add rownames and colnames
        res <- data.frame(Row = rep(1:nrow(rMinv), nrow(rMinv)),
                          Column = rep(1:nrow(rMinv), each = nrow(rMinv)),
                          coeff = as.numeric(rMinv),
                          lower = as.logical(lower.tri(rMinv, diag = TRUE)))
        rm(rMinv)
                        
      
    
        # only use lower triangle
        res <- res[res$lower == TRUE, c("Row", "Column", "coeff")]
          
        if (sorting=="ASReml"){    
          res <-  res[order( res$Row,  res$Column), ] 
        }
        if (sorting=="WOMBAT"){
          res <- res[, c(2,1,3)]
          res <-  res[order(res$Column,  res$Row), ]  
        }
        res <- res[res$coeff != 0, ]
        
        # write to given file
        if (!is.null(file)){
          write.table(res, file, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
          rm(x)
        } else {
          attr(res, "rowNames") <- rownames(x)
          rm(x)
          return(res)
        }

}
