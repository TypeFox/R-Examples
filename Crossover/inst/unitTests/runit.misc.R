test.design.functions <- function () {
  s <- 3
  p <- 4
  v <- 3
  
  v.rep <- rep((s*p) %/% v, v) + c(rep(1, (s*p) %% v), rep(0, v-((s*p) %% v)))
  design <- matrix(sample(rep(1:v, v.rep)), p, s)
  
  rcDesign <- rcd(design, v, model=1)
  # JRW, p 2650, first equation on that page, whithout number
  Ar <- infMatrix(rcDesign, v, model=1)
  Xr <- rcdMatrix(rcDesign, v, model=1)
  # JRW, p 2650, second equation on that page, number 11
  Ar2 <- t(Xr) %*% (diag(s*p)-getPZ(s,p)) %*% Xr
  checkTrue(max(abs(Ar-Ar2))<0.00001)
  
  Csub <- contrMat(n=rep(1, v), type="Tukey")
  class(Csub) <- "matrix" #TODO Package matrix can be improved here (IMO)!
  C <- as.matrix(bdiag(Csub,Csub))
  H <- linkMatrix(model=1, v)
  var1 <- sum(diag(C %*% ginv(t(H) %*% Ar %*% H) %*% t(C)))
  
  gco <- general.carryover(design, model=1)
  var2 <- sum(gco$Var.trt.pair[lower.tri(gco$Var.trt.pair)]) + sum(gco$Var.car.pair[lower.tri(gco$Var.car.pair)])
  
  # checkTrue(abs(var1-var2)<0.00001) This is often not true due to the fact that s*p are not much bigger than v.
}

test.brute.force.compare.approaches <- function() {
  if (!"extended" %in% strsplit(Sys.getenv("CROSSOVER_UNIT_TESTS"),",")[[1]]) {
    cat("Skipping design tests for comparing approaches.\n")
    return()
  }
  
  for (model in c(1:6,8)) { #TODO model 7 differs for general.carryover and search algorithm
    s <- 10 #TODO Also make s, p and v random.
    p <- 6
    v <- 3  
    v.rep <- rep((s*p) %/% v, v) + c(rep(1, (s*p) %% v), rep(0, v-((s*p) %% v)))
    
    differences <- 0
    
    for (i in 1:10) {
      cat("\nRun ",i,":\n")
      design <- matrix(sample(rep(1:v, v.rep)), p, s)
      result <- try(Crossover:::compareApproaches(design, models2check=model, stop.on.diff=TRUE))
      if ("try-error" %in% class(result)) {
        differences <- differences + 1
      }
    }
    if(differences>5) {
      stop(paste("general.carryover and search algorithm differ in ",differences, " out of 10 cases for model ",model,".", sep=""))
    }
  }  
}

test.strangeDesignInputs <- function() {
  s <- 4 # number of sequences
  p <- 4 # number of periods
  v <- 4 # number of treatments
  
  D <- rbind(c("A","B","C","D"),
             c("B","C","D","A"),
             c("C","D","A","B"),
             c("D","A","B","C"))
  
  D <- matrix(as.numeric(as.factor(D)), dim(D)[1])  
  
  myInv <- ginv(rcd(D, v, model=1))
  
}

test.getEff <- function() {
  checkTrue(all(abs(getEff(getDesign("pidgeon1"))-c(0.712893817102914, 0.712893817102914, 0.715489015631601, 0.712893817102914, 
                                                    0.712893817102914, 0.87962962962963, 0.106076388888889, 0.201785483035483, 
                                                    0.87962962962963))<0.00001))
}