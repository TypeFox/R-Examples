`summary.segmented` <-
function(object, short=FALSE, var.diff=FALSE, ...){
    if(is.null(object$psi)) object<-object[[length(object)]]
#i seguenti per calcolare aa,bb,cc funzionano per lm e glm, da verificare con arima....
#    nome<-rownames(object$psi)
#    nome<-as.character(parse("",text=nome))
#    aa<-grep("U",names(coef(object)[!is.na(coef(object))]))
#    bb<-unlist(sapply(nome,function(x){grep(x,names(coef(object)[!is.na(coef(object))]))},simplify=FALSE,USE.NAMES=FALSE))
#    cc<-intersect(aa,bb) #indices of diff-slope parameters
#    iV<- -grep("psi.",names(coef(object)[!is.na(coef(object))]))#indices of all but the Vs
    if(var.diff && length(object$nameUV$Z)>1) {
      var.diff<-FALSE
      warning("var.diff set to FALSE with multiple segmented variables", call.=FALSE)
      }
    nomiU<-object$nameUV[[1]]
    nomiV<-object$nameUV[[2]]
    idU<-match(nomiU,names(coef(object)[!is.na(coef(object))]))
    idV<-match(nomiV,names(coef(object)[!is.na(coef(object))]))
    beta.c<- coef(object)[nomiU]
    #per metodo default..
    if( !inherits(object, "segmented")){
      summ <- c(summary(object, ...), object["psi"])
      summ[c("it","epsilon")]<-object[c("it","epsilon")]
      coeff<-coef(object)
      v<-try(vcov(object), silent=TRUE)
      if(class(v)!="try-error"){
          v<-sqrt(diag(v))
          summ$gap<-cbind(coeff[idV]*beta.c,abs(v[idV]*beta.c),coeff[idV]/v[idV])
          colnames(summ$gap)<-c("Est.","SE","t value")
          rownames(summ$gap)<-nomiU
      } else {
          summ$gap<-cbind(coeff[idV]*beta.c,NA,NA)
          colnames(summ$gap)<-c("Est.","SE","t value")
          rownames(summ$gap)<-nomiU
      }
    return(summ)
      }
    if("lm"%in%class(object) && !"glm"%in%class(object)){
    #if(!inherits(object, "glm")){
        summ <- c(summary.lm(object, ...), object["psi"])
        summ$Ttable<-summ$coefficients
        if(var.diff){
            sigma2.new<-tapply(object$residuals, object$id.group, function(xx){sum(xx^2)}) 
            summ$df.new<-tapply(object$residuals, object$id.group, function(xx){(length(xx)-length(object$coef))})
            summ$sigma.new<-sqrt(sigma2.new/summ$df.new)
            #modifica gli SE
            Qr <- object$qr
            p <- object$rank
            p1 <- 1L:p
            inv.XtX <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
            X <- qr.X(Qr,FALSE)
            attr(X, "assign") <- NULL
            sigma.i<-rowSums(model.matrix(~0+factor(object$id.group))%*%diag(summ$sigma.new))
            var.b<-inv.XtX%*%crossprod(X*sigma.i)%*%inv.XtX
            dimnames(var.b)<-dimnames(summ$cov.unscaled)
            summ$cov.var.diff<-var.b
            summ$Ttable[,2]<-sqrt(diag(var.b))
            summ$Ttable[,3]<-summ$Ttable[,1]/summ$Ttable[,2]
            summ$Ttable[,4]<- 2 * pnorm(abs(summ$Ttable[,3]), lower.tail = FALSE)
            dimnames(summ$Ttable) <- list(names(object$coefficients)[Qr$pivot[p1]],
                  c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
          }
        coeff<-summ$Ttable[,1]#summ$coefficients[,1]
        v<-summ$Ttable[,2] #summ$coefficients[,2]
        summ$gap<-cbind(coeff[idV]*beta.c,abs(v[idV]*beta.c),coeff[idV]/v[idV])
        summ$Ttable[idU,4]<-NA
        summ$Ttable<-summ$Ttable[-idV,] 
        #dimnames(summ$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        colnames(summ$gap)<-c("Est.","SE","t value")
        rownames(summ$gap)<-nomiU
        summ[c("it","epsilon")]<-object[c("it","epsilon")]
        summ$var.diff<-var.diff
        summ$short<-short
        class(summ) <- c("summary.segmented", "summary.lm")
        return(summ)
        }
    #if("glm"%in%class(object)){
    if(inherits(object, "glm")){
        summ <- c(summary.glm(object, ...), object["psi"])
        summ$Ttable<-summ$coefficients[-idV,]
        summ$Ttable[idU,4]<-NA
        coeff<-summ$coefficients[,1]
        v<-summ$coefficients[,2]
        summ$gap<-cbind(coeff[idV]*beta.c,abs(v[idV]*beta.c),coeff[idV]/v[idV])
        #dimnames(summ$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        colnames(summ$gap)<-c("Est.","SE","t value")
        rownames(summ$gap)<-nomiU
        summ[c("it","epsilon")]<-object[c("it","epsilon")]
        summ$short<-short
        class(summ) <- c("summary.segmented", "summary.glm")
        return(summ)}
    if("Arima"%in%class(object)){
        #da controllare
        coeff<-object$coef
        v<-sqrt(diag(object$var.coef))
        Ttable<-cbind(coeff[-idV],v[-idV],coeff[-idV]/v[-idV])
        object$gap<-cbind(coeff[idV]*beta.c,v[idV]*beta.c,coeff[idV]/v[idV])
        #dimnames(object$gap)<-list(rep("",nrow(object$psi)),c("Est.","SE","t value"))
        colnames(object$gap)<-c("Est.","SE","t value")
        rownames(object$gap)<-nomiU
        colnames(Ttable)<-c("Estimate","Std. Error","t value")
        object$Ttable<-Ttable
        object$short<-short
        summ<-object 
        class(summ) <- c("summary.segmented", "summary.Arima")
        return(summ)}
}

