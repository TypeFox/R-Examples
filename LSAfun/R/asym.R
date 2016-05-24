## Asymmetrical similarity measures 

#' @export

asym <- function(x,y,method,t=0,tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    if(breakdown==TRUE){
      
      x <- breakdown(x)      
      y <- breakdown(y)
      
    }
    
    
    if(x %in% rownames(tvectors) && y %in% rownames(tvectors)){
      
      v <- tvectors[x,]
      w <- tvectors[y,]
      
          
      
      if(method == "weedsprec"){
        
        ## compute active dimensions
        
        if((any(v<0) | any(w<0)) && t < 0){
          
          warning("vector for x or y contains negative values, results will be flawed")
          
        }
        
        if(all(v <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        if(all(w <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        active_v <- which(v > t)
        active_w <- which(w > t)
        
        ## identify simultaneously active dimensions
        active_vw <- intersect(active_v,active_w)
        
        ## compute weedsprec
        weedsprec <- sum(v[active_vw])/sum(v[active_v])
        return(weedsprec)
        
      }
      
      
      if(method == "cosweeds"){
        
        ## compute active dimensions
        
        if((any(v<0) | any(w<0)) && t < 0){
          
          warning("vector for x or y contains negative values, results will be flawed")
          
        }
        
        if(all(v <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        if(all(w <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        
        active_v <- which(v > t)
        active_w <- which(w > t)
        
        ## identify simultaneously active dimensions
        active_vw <- intersect(active_v,active_w)
        
        ## compute weedsprec and cosine
        weedsprec <- sum(v[active_vw])/sum(v[active_v])
        
        cos       <- sum(v[active_vw]*w[active_vw])/
          ( sqrt( sum((v[active_v])^2) ) * sqrt( sum((w[active_w])^2) ) )  
        
        ## compute cosweeds
        
        cosweeds <- sqrt(weedsprec*cos)
        return(cosweeds)
        
      }
      
      
      if(method == "clarkede"){
        
        ## compute active dimensions
        
        if((any(v<0) | any(w<0)) && t < 0){
          
          warning("vector for x or y contains negative values, results will be flawed")
          
        }
        
        if(all(v <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        if(all(w <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        
        active_v <- which(v > t)
        active_w <- which(w > t)
        
        ## identify simultaneously active dimensions
        active_vw <- intersect(active_v,active_w)
        
        ## compute clarkede
        av   <- v[active_vw]
        aw   <- w[active_vw]
        list <- data.frame(av,aw)
        mins <- apply(list,1,min)
        
        clarkede <- sum(mins)/sum(v[active_v])
        return(clarkede)
        
      }
      
      
      if(method == "invcl"){
        
        ## compute active dimensions
        
        if((any(v<0) | any(w<0)) && t < 0){
          
          warning("vector for x or y contains negative values, results will be flawed")
          
        }
        
        
        if(all(v <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        if(all(w <= t)){
          
          stop("no entry of the vector for x is active (all values < t or = t)")
        }
        
        
        active_v <- which(v > t)
        active_w <- which(w > t)
        
        ## identify simultaneously active dimensions
        active_vw <- intersect(active_v,active_w)
        
        ## compute clarkede
        av   <- v[active_vw]
        aw   <- w[active_vw]
        list <- data.frame(av,aw)
        mins <- apply(list,1,min)
        
        clarkede <- sum(mins)/sum(v[active_v])
        
        ## compute invcl
        invcl <- sqrt(clarkede*(1-clarkede))
        return(invcl)
        
      }
      
      
      if(method=="kintsch"){
        
        LV <- sqrt(sum(v^2))
        LW <- sqrt(sum(w^2))
        
        cos <- as.numeric(cosine(v,w))
        
        V_W <- LW*cos
        W_V <- LV*cos
        
        SIM_VW <- V_W/max(LV,LW)
        return(SIM_VW)
        
        
      }
      
      
      
      
    }else{
      
      warning("x and y must be words in rownames(tvectors)")
      return(NA)
      
    }
    
  }else{stop("tvectors must be a matrix!")}
  
  
  
}