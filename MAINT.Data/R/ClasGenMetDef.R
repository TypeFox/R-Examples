setClass("IData",representation(MidP="data.frame",LogR="data.frame",ObsNames="character",VarNames="character",NObs="numeric",NIVar="numeric"))
setClass("IdtE",representation(ModelNames="character",ModelType="character",ModelConfig="numeric",NIVar="numeric",SelCrit="character",
				logLiks="numeric",AICs="numeric",BICs="numeric",BestModel="numeric"),contains="VIRTUAL")
setClass("IdtSngDE",contains=c("IdtE","VIRTUAL"))
setClass("IdtMxE",representation(grouping="factor"),contains=c("IdtE","VIRTUAL"))
setClass("IdtSngNDE",representation(mleNmuE="numeric",mleNmuEse="numeric",Configurations="list"),contains="IdtSngDE")
setClass("IdtMxNDE",representation(Hmcdt="logical",mleNmuE="matrix",mleNmuEse="matrix",Configurations="list"),contains="IdtMxE")
setClassUnion("IdtNDE",c("IdtSngNDE","IdtMxNDE"))
setClass("LRTest",representation(QuiSq="numeric",df="numeric",pvalue="numeric",H0logLik="numeric",H1logLik="numeric"))
setClass("ConfTests",representation(TestRes="list",RestModels="character",FullModels="character"))
setClass("IdtMANOVA",representation(NIVar="numeric",grouping="factor",H0res="IdtSngDE",H1res="IdtMxE"),contains="LRTest")
setClass("IdtClMANOVA",contains="IdtMANOVA")
setClass("IdtHetNMANOVA",contains="IdtMANOVA")
setClass("Idtlda",representation(prior="numeric",means="matrix",scaling="matrix",N="numeric"))
setClass("Idtqda",representation(prior="numeric",means="matrix",scaling="array",ldet="numeric",lev="character"))

setGeneric("nrow")
setGeneric("ncol")
setGeneric("summary",signature="object")
setGeneric("head",package="utils",signature="x")
setGeneric("tail",package="utils",signature="x")
setGeneric("coef",package="stats",signature="object")
setGeneric("predict",package="stats",signature="object")

setGeneric("stdEr",package="miscTools",signature="x")
setGeneric("lda",package="MASS",signature="x")
setGeneric("qda",package="MASS",signature="x")

setGeneric("mle",function(Idt,Model="Normal",Config=1:5,SelCrit=c("AIC","BIC"))  standardGeneric("mle"))
setGeneric("MANOVA",function(Idt, grouping, Model="Normal", Config=1:5, SelCrit=c("AIC","BIC"), Mxt=c("Hom","Het"), tol=1.0e-4)  standardGeneric("MANOVA"))
setGeneric("BestModel",function(Idt,SelCrit=c("IdtCrt","AIC","BIC"))  standardGeneric("BestModel"))
setGeneric("testMod",function(Idt,RestMod=1:length(Idt@ModelNames),FullMod="Next")  standardGeneric("testMod"))
setGeneric("H1res",  function(object) standardGeneric("H1res"))
setGeneric("H0res",  function(object) standardGeneric("H0res"))



