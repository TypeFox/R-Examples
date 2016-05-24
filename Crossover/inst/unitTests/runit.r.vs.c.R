test.r.vs.c <- function () {
  s <- 3
  p <- 4
  v <- 3
  for (model in 1:8) {      
      v.rep <- rep((s*p) %/% v, v) + c(rep(1, (s*p) %% v), rep(0, v-((s*p) %% v)))
      design <- matrix(sample(rep(1:v, v.rep)), p, s)
      
      rcDesign <- rcd(design, v, model=model)
      rcDesign_R <- rcd_R(design, v, model=model)
      checkTrue(max(abs(rcDesign-rcDesign_R))<0.00001)
      
      options(warn=2)
      checkException(rcdMatrix(design, v, model=model))
      checkException(rcdMatrix_R(design, v, model=model))
      options(warn=1)
      
      rcDesignMat <- rcdMatrix(rcDesign, v, model=model)
      rcDesignMat_R <- rcdMatrix_R(rcDesign, v, model=model)
      checkTrue(max(abs(rcDesignMat-rcDesignMat_R))<0.00001)
      
      # TODO Implement infMatrix for model 3 and 7:
      if (!model %in% c(3,7)) {
        Ar <- infMatrix(rcDesign, v, model=model)
      } else {
        Ar <- infMatrix_R(rcDesign, v, model=model, method=1)  
      }
      Ar_R1 <- infMatrix_R(rcDesign, v, model=model, method=1)  
      Ar_R2 <- infMatrix_R(rcDesign, v, model=model, method=2)  
      checkTrue(max(abs(Ar-Ar_R1))<0.00001)
      checkTrue(max(abs(Ar_R2-Ar_R1))<0.00001)
      
      d <- structure(c(2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 2, 
                       2, 1, 1, 2, 1, 1, 2), .Dim = c(4L, 6L))
      checkTrue(estimable_R(d, v=2, model=2)==estimable(d, v=2, model=2))
  }
}