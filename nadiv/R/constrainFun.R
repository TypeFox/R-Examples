constrainFun <- function(parameter.val, full, fm2, comp, G, mit = 600){
  row <- which(fm2$Gamma == comp)
  fm2[row, 2:3] <- c(parameter.val, "F")
  if(G) full$G.param <- fm2 else full$R.param <- fm2
  con.mod <- asreml::update.asreml(object = full, maxiter = mit, trace = FALSE)
  cnt <- 0
  while(!con.mod$converge & cnt <= 5){
     con.mod <- asreml::update.asreml(con.mod)
     cnt <- cnt + 1
  }
  cnt <- 0
  if(con.mod$converge){
     pcc.out <- pcc(con.mod, silent = TRUE)
     while(!pcc.out & cnt <= 5){
        con.mod <- asreml::update.asreml(con.mod, maxiter = mit)
        if(con.mod$converge) pcc.out <- pcc(con.mod, silent = TRUE)
        cnt <- cnt + 1
     }
     con.mod$converge <- pcc.out
  }
  if(con.mod$converge) return(LRTest(full$loglik, con.mod$loglik)$lambda) else return(NA)
}

