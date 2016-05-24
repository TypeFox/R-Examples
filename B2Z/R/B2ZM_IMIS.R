B2ZM_IMIS <- function(data, priorBeta, priorQ, priorG, v, S, 
                      tauN.sh, tauN.sc, tauF.sh, tauF.sc,  
                      VN, VF, indep.model = FALSE, cred = 95,
                      N0 = 6000, B = 600, M = 3000, it.max = 16, 
                      figures = list(save = FALSE, 
                      type =c("ps", "eps","pdf", "png", "jpg"))) {

   if(!indep.model){
      tauN.sh = tauN.sc = tauF.sh = tauF.sc = 0
   }
   else{
      S <- matrix(0, 2, 2)
      v <- 0
   }

   ans <- B2ZM(data = data, priorBeta = priorBeta, priorQ = priorQ, priorG = priorG,
               v = v, S = S, tauN.sh = tauN.sh, tauN.sc = tauN.sc, tauF.sh = tauF.sh, 
               tauF.sc = tauF.sc, VN = VN, VF = VF, indep.model = indep.model, 
               cred = 95,  sampler = "IMIS", imis.control = list(N0 = N0, 
               B = B, M = M, it.max = it.max),  figures = figures)

   return(ans)
   }