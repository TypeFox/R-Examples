require(lmerTest)

load(system.file("testdata","NSstdData.RData", package="lmerTest"))


Data.til.R$Consumer <- as.factor(Data.til.R$Consumer)
Data.til.R$Age <- as.factor(Data.til.R$Age)
Data.til.R$PresPosition <- as.factor(Data.til.R$PresPosition)



model <- lmer(Balanced ~ Type*Exposure +(1|PresPosition), data=Data.til.R)
c <- step(model, reduce.random=FALSE)


stopifnot(dim(c$lsmeans.table)[1] > 0)

