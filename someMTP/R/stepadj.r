
step.adj<-function(object, MC=1000, scope = NULL, scale = 0, direction = c("both", "backward", "forward"), trace = 0,keep = NULL,steps = 1000, k = 2)
{
family=object$family
if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
object$family=family
TEST=ifelse((object$family[1]=="gaussian")|(is(object,"lm")),"F","Ch")
direction <- match.arg(direction)

p=rep(NA,MC)

#perform null model and compute p-value
formula0=object$formula
formula0[[3]]=1
mod0=glm(formula=formula0,family=object$family, control =object$control,data=object$data)
mod.step=step(object,scope=scope, scale=scale,direction=direction ,trace=trace,keep=keep,steps=steps, k=k )
anova.sw=anova(mod0,mod.step,test=TEST)
p[MC]=anova.sw[[6]][2]

#create redo function that sample the response, re-fit the model and step() select.
redo<-function(mod,mod0,TEST,scope, scale,direction,keep,steps, k){
	mod$data[,1]=sample(mod$data[,1])
	mod=update(mod,data=mod$data)
	mod.step=step(mod,scope=scope, scale=scale,direction=direction ,trace=0,keep=keep,steps=steps, k=k )
	anova.swr=anova(mod0,mod.step,test=TEST)[[5+(TEST=="F")]][2]
	p=min(anova.swr,1,na.rm =T)
}
p[-MC]=replicate(MC-1,redo(object,mod0,TEST,scope, scale,direction,keep,steps, k))


#compute the corrected p-value using mean(p<=p[MC])
attributes(anova.sw)$heading[1]=paste("Analysis of Deviance Table on Selected Model \nwith Multiplicity Correction (weak FWE-control)\nbased on",MC,"Random permutations\n")
anova.sw[[" sw.adj P"]]=c(NA, mean(p<=p[MC]))

return(anova.sw=anova.sw)
}
