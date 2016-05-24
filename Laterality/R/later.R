later <- function (data, catch="Food", indiv="Indiv", hand="Hand", RightHand = "R", LeftHand = "L", savetable = FALSE, file="HIz.csv")
{
  for (i in 1:nlevels(data[[catch]])) {
      HIobs<-data[data[[catch]]==levels(data[[catch]])[i],]
      Tab<-table(HIobs[[indiv]],HIobs[[hand]])
      Dat<-as.data.frame.matrix(Tab) #contingency table
      
      ifelse (is.null(Dat[[RightHand]]) == TRUE, HI<-(-Dat[[LeftHand]])/Dat[[LeftHand]], ifelse (is.null(Dat[[LeftHand]]) == TRUE, HI<-Dat[[RightHand]]/Dat[[RightHand]], HI<-(Dat[[RightHand]]-Dat[[LeftHand]])/(Dat[[RightHand]]+Dat[[LeftHand]]))) #Handedness index      
      Dat$HI<-HI #add HI to the table
      
      ifelse (is.null(Dat[[RightHand]]) == TRUE, z<-(-Dat[[LeftHand]]/2)/sqrt(Dat[[LeftHand]]/4), ifelse (is.null(Dat[[LeftHand]]) == TRUE, z<-(Dat[[RightHand]]-Dat[[RightHand]]/2)/sqrt(Dat[[RightHand]]/4), z<-(Dat[[RightHand]]-(Dat[[LeftHand]]+Dat[[RightHand]])/2)/sqrt((Dat[[LeftHand]]+Dat[[RightHand]])/4)))  #z-score
      Dat$z<-z #add z to the table
	    
      p.val<-2*(pnorm(abs(z), lower.tail=FALSE))
      Dat$p.val<-p.val #add the p.val to the table
      
      reshand <- "A"
      reshand <- ifelse (z>=1.96, "R", reshand)
      reshand <- ifelse (z<=-1.96, "L", reshand)
      Dat$Hand<-reshand #add Hand to the table
      
      signif <- "."
      signif <- ifelse (p.val<=0.05, "*", signif)
      signif <- ifelse (p.val<=0.01, "**", signif)
      Dat$Signif<-signif #add Signif to the table
      
      if("HIz" %in% ls() == FALSE) {
        HIz<-matrix(,nlevels(data[[indiv]]),0)
      } else {
        }
      
      HIz<-cbind(HIz,Dat)
  }
  
  HImatrix<-as.matrix(HIz)
  HIzar<-array(HImatrix,dim=c(nlevels(data[[indiv]]),ncol(Dat),nlevels(data[[catch]])))
  dimnames(HIzar) <- list(c(levels(data[[indiv]])), c(levels(data[[hand]]), "HI", "z", "p.val", "Hand", "Signif"), c(levels(data[[catch]])))
  
  if (savetable == "csv") {write.csv(HIzar[,,], file = file)} else {}
  if (savetable == "csv2") {write.csv2(HIzar[,,], file = file)} else {}
  
  message("-----------------------------------------")       
  message("-1 </= HI </= 1")
  message("If HI > 0 -> right hand preference")
  message("If HI < 0 -> left hand preference")
  message("-----------------------------------------")       
  message("for alpha = 0.05")
  message("If z </= -1.96 -> left-handed (L)")
  message("If z >/= 1.96 -> right-handed (R)")
  message("If -1.96 < z < 1.96 -> ambiguously handed (A)")  
  message("-----------------------------------------")
  message("if p.val > 0.05 -> .")
  message("if p.val <= 0.05 -> *")
  message("if p.val <= 0.01 -> **")
  message("-----------------------------------------")
         
  HIzar
}

