summary.FindIt<-function(object,...){
    ## x <- summary.glm(object, dispersion = dispersion,
    ## correlation = correlation, symbolic.cor = symbolic.cor, ...)
    treat.type <- object$treat.type
    type <- object$type
    main <- object$main
    not.zero<-object$coefs.orig!=0
    
    if(type=="binary"){
        coef.table<-cbind(object$names.out,object$coefs.orig/2)[not.zero,]
        coef.print<-cbind(object$names.out,signif(object$coefs.orig/2,3))[not.zero,]
    }
    if(type=="continuous"){
        coef.table<-cbind(object$names.out,object$coefs.orig)[not.zero,]
        coef.print<-cbind(object$names.out,signif(object$coefs.orig,3))[not.zero,]
    }
    
    rownames(coef.print)<-rownames(coef.table)<-object$names.out[not.zero]	
    coef.print[coef.print=="0"]<-"0.000"	
    colnames(coef.print)<-c("Coefficient","Estimate")

    null.gcv<-var(object$y)/(length(object$y)-1)*length(object$y)
    model.gcv<-(object$GCV)*length(object$y)

    if(main){
        model.main  <- object$model.main
    }
    if(main & treat.type=="single"){
        model.int   <- object$model.int
    }
    
    model.treat <- object$model.treat
    
    cat("\nCall:\n")
    cat(" Treatment Model: ")
    print( model.treat)
    
    if(main){
        cat(" Main Model : ")
        print(model.main)
    }
    
    if(main & treat.type=="single"){
        cat(" Interaction Covariates: ")
        print(model.int)
    }
    

    cat(" Treatment type: ")
    print(treat.type)

    cat(" Outcome type: ")
    print(type)
    
    
    cat("\nATE:\n")
    print(object$ATE)

    
    cat("\nCoefficients:\n")

    print(noquote(coef.print[,-1]))  
    
    cat("\n---------")
    
    cat("\nModel Fit Statistics:\n")
    cat(c("GCV:\n"))
    cat(c("  Null: ",round(null.gcv,3)))
    cat(c("  Model: ",round(model.gcv,3),"\n"))
    cat(c("Percent Misclassified:\n"))
    cat(c("   Null: ",round(min(c(mean(sign(object$y)==-1),mean(sign(object$y)==1))),2)))
    cat(c("   Model: ",round(mean(sign(object$y)!=sign(object$fit)),2),"\n"))
    cat(c("   Percent Improvement, vs. NULL: ",round(100-100*mean(sign(object$y)!=sign(object$fit))/min(c(mean(sign(object$y)==-1),mean(sign(object$y)==1))),2),"% \n"))
    cat(c("Percent Outside Margin:\n ",round(mean((object$y^2- object$y*object$fit)<=0 )*100,3),"%, n =",sum((1- object$y*object$fit)<=0 ), "\n")) 	
    out<-list("coefficients"=noquote(coef.print[,-1]),
              "GCV"=c(null.gcv,model.gcv),
              "misclass"=c(min(c(mean(sign(object$y)==-1),mean(object$y==1))),mean(sign(object$y)!=sign(object$fit))))
    invisible(out)  
}
