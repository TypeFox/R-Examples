B2ZM_SIR <- function(data, priorBeta, priorQ, priorG, 
                     v, S, tauN.sh, tauN.sc, tauF.sh, tauF.sc,  VN,
                     VF, indep.model = FALSE, cred = 95,  m = 10000, 
                     figures = list(save = FALSE, 
                     type =c("ps", "eps","pdf", "png", "jpg"))) {


   if(!indep.model){
      tauN.sh = tauN.sc = tauF.sh = tauF.sc = 0
   }
   else{
      S <- matrix(0,2,2)
      v <- 0
   }

   ans <- B2ZM(data = data, priorBeta = priorBeta, priorQ = priorQ,  
               priorG = priorG, v = v, S = S, tauN.sh = tauN.sh, 
               tauN.sc = tauN.sc, tauF.sh = tauF.sh, tauF.sc = tauF.sc, 
               VN = VN, VF = VF, indep.model = indep.model, cred = 95,  
               sampler = "SIR",  sir.control = list(m = m), figures = figures)

   return(ans)
   }