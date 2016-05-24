
mse<-function(repstat, design){
       v<-attr(repstat,"var")
       center<-attr(v,"means")
       if ((length(v)!=length(center)^2) && (length(v)==length(center))){
         attr(repstat,"var")<-vcov(repstat)+(center-coef(repstat))^2*sum(design$rscales)*design$scale
       } else {
         attr(repstat,"var")<-as.matrix(vcov(repstat)+outer((center-coef(repstat)))*sum(design$rscales)*design$scale)
       }
       repstat
       }

