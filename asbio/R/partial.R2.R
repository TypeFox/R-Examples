partial.R2<-function(nested.lm, ref.lm){
a <- anova(nested.lm)
b <- anova(ref.lm)
length.ref <- length(attributes(ref.lm$terms)$"dataClasses")
length.nested <- length(attributes(nested.lm$terms)$"dataClasses")
if(length.nested > length.ref) stop("Specify nested model first in arguements")
if(length.ref - length.nested > 1) stop("Reference and nested model should only differ with repsect to the presence/absence of one predictor")
SSE.wo <- tail(a$"Sum Sq", 1)
SSE.with <- tail(b$"Sum Sq", 1)
P.R2<-(SSE.wo-SSE.with)/SSE.wo
P.R2
}


