## SIM in context - function following Kintsch (2014)

#' @export

conSIM <- function(x,y,z,c,tvectors=tvectors,breakdown=FALSE){
  
  if(class(tvectors) == "matrix"){
    
    if(breakdown==TRUE){
      
      x <- breakdown(x) 
      y <- breakdown(y)
      z <- breakdown(z)
      c <- breakdown(c)
      
    }
    
    
    if(x %in% rownames(tvectors) && y %in% rownames(tvectors)){
      
      A <- tvectors[x,]
      B <- tvectors[y,]
      C <- tvectors[z,]
      X <- tvectors[c,]
      
      LA <- sqrt(sum(A^2))
      LB <- sqrt(sum(B^2))
      LC <- sqrt(sum(C^2))
      LX <- sqrt(sum(X^2))
      
      
      ## aux function
      
      minisim <- function(LA,LB,cos){
        
        A_B <- LB*cos
        B_A <- LA*cos
        
        SIM_AB <- A_B/max(LA,LB)
        SIM_BA <- B_A/max(LA,LB)
        
        result <- list(SIM_xy=SIM_AB,SIM_yx=SIM_BA)
        return(result)
        
      }
      
      
      ## first
      
      cosAB <- as.numeric(cosine(A,B))
      
      CX  <- C+X
      LCX <- sqrt(sum(CX^2)) 
      
      cosBCX <- as.numeric(cosine(B,CX))
      LBCX   <- LCX*cosBCX
      
      LB_star <- LB - LBCX
      
      SIM_A_BCX_star <- minisim(LA,LB_star,cosAB)$SIM_xy
      
      
      
      
      ## second
      
      cosAC <- as.numeric(cosine(A,C))
      
      BX  <- B+X
      LBX <- sqrt(sum(BX^2)) 
      
      cosCBX <- as.numeric(cosine(C,BX))
      LCBX   <- LBX*cosCBX
      
      LC_star <- LC - LCBX
      
      SIM_A_CBX_star <- minisim(LA,LC_star,cosAC)$SIM_xy
  
      
      result <- list(SIM_XY_zc=SIM_A_BCX_star,SIM_XZ_yc=SIM_A_CBX_star)
      
      return(result)
      
      
    }else{
      
      warning("x and y must be words in rownames(tvectors)")
      return(NA)
      
    }
    
  }else{stop("tvectors must be a matrix!")}
  
}