setOldClass("topmod") #,representation(addmodel="function", lik="function", bool="function", ncount="function", nbmodels="function", nregs="function", betas_raw="function", betas2_raw="function", kvec_raw="function", bool_binary="function", betas="function", betas2="function", fixed_vector="function"))
setOldClass("bma") #,representation(info="list",arguments="list",topmod="topmod",start.pos="integer",gprior.info="list",mprior.info="list",X.data="data.frame",reg.names="character",bms.call="call"))
setOldClass("zlm") #,representation(coefficients="numeric",residuals="numeric",rank="numeric",fitted.values="numeric",df.residual="numeric",call="call",terms="formula",model="data.frame",coef2moments="numeric",marg.lik="numeric",gprior.info="list"))


#setMethod("[","bma",.index.bma)
#setMethod("c","bma",c.bma)
#setMethod("coef","bma",coef.bma)
#setMethod("density","bma",density.bma)
#setMethod("deviance","bma",deviance.bma)
#setMethod("image","bma",image.bma)
#setMethod("model.frame","bma",model.frame.bma)
#setMethod("plot","bma",plot.bma)
#setMethod("predict","bma",predict.bma)
#setMethod("print","bma",print.bma)
#setMethod("summary","bma",info.bma)
#setMethod("variable.names","bma",variable.names.bma)

#zlm Methods
#setMethod("density","zlm",density.zlm)
#setMethod("deviance","zlm",deviance.zlm)
#setMethod("logLik","zlm",logLik.zlm)
#setMethod("predict","zlm",predict.zlm)
#setMethod("summary","zlm",summary.zlm)
#setMethod("variable.names","zlm",variable.names.zlm)
#setMethod("vcov","zlm",vcov.zlm)

#topmod methods
#setMethod("print","topmod",print.topmod)
#setMethod("[","topmod",.index.topmod)