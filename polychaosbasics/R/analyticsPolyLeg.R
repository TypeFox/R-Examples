###################################################
# analyticsPolyLeg
# The input dataset is built by the function 'randomLHS'
# from package 'lhs', then response is calculated
# from standard models 

analyticsPolyLeg <- function(nlhs,degree,
                         model.fun)

{
  # NOTE: result depends on the random seed
  if (degree <= 1) {
      stop("The degree should be greater than 1")
    }
  
  # Initialisations
  switch (model.fun,
         ishigami = {
           nvx <- 3
         },
          sobol= {
           nvx <- 4
         },
          polymodel = {
              # removed from this version
               stop("Invalid argument 'model.fun'. Valid values are 'ishigami', 'sobol'")
            nvx <- 3
         },
          stop("Invalid argument 'model.fun'. Valid values are 'ishigami', 'sobol', 'polymodel'")
          ) # end switch
  
 Y <- NULL 
# Construction du plan lhd a valeurs dans [0,1]
 lhd0 <-randomLHS(nlhs,nvx)
  
 if (model.fun == "sobol") {
 # Rajout des  variables fixes
   vfix <- matrix(0.5, nrow=nlhs, ncol=4)
   planInt<-cbind(lhd0,vfix)
   Y<-sobolModel(planInt)
 }
 
# Passage au rang: [1,nlhs]
 lhd1<-apply(lhd0,2,rank,ties.method = "random")


#  Calcul du model (des Yi)
 switch (model.fun,
         ishigami = {
          # Lhd "naturel" = dans les "vrais" bornes
           binf <- matrix(-pi,nrow=nvx,ncol=1)
           bsup<-matrix(pi,nrow=nvx,ncol=1)
           lhdnat<-calibLhd(lhd1,binf,bsup)
           # Values of fixed model factors 
           para<-c(7,2,0.1,4)
           Y<-IshigamiModel(lhdnat, para)},
         sobol= {
           # Y already calculated
         },
         polymodel = {
           Y<-polyModel(lhd0)
         }
         ) #end switch
 
#############################################
#  Construction des monomes
#############################################
 plan2<-Structure(nvx,degree)
  dimnames(plan2) <- list(c("0", labelmono(plan2)),
                          paste("Input", 1:nvx))

# Calibration du lhd de base sur la plage [-1,1]
# pour polynomes de Legendre
 binfl<-rep(-1,nvx) ; bsupl<-rep(1,nvx)
# Appel calibLhd()
 lhdc<-calibLhd(lhd1,binfl,bsupl)
#############################################
#  Construction de la matrice du modele 
#############################################
 XM<-modLeg(lhdc,degree,plan2); XM<-as.matrix(XM)
#############################################
# return
#############################################
XMY <- cbind(XM,Y)

   retour <- new("PCEpoly", .Data=XMY, design=plan2, nvx=nvx,
                 call=match.call())
 return(retour)
}

