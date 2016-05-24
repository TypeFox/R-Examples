library(gRbase)
data(HairEyeColor)
gm.hec <- as.gmData(HairEyeColor)

varNames(gm.hec)
print(gm.hec)

m1i <- hllm(~Hair+Eye+Sex, gm.hec) 
m2i <- hllm( ~Hair*Eye+Sex, gm.hec) 

m1f <- fit(m1i); class(m1f)
m2f <- fit(m2i); class(m2f)

m1i <- hllm( ~. , gmData=gm.hec) 
m2i <- hllm( ~., gm.hec) 

m1s <- hllm( ~.. , gmData=gm.hec) 
m2s <- hllm( ~.. , gmData=gm.hec) 


### Test of several models
mod.list <- c(~., ~.., ~..., ~...., ~.^1, ~.^2, ~.^3, ~.^4, ~Hair+Eye+Sex, ~Hair*Eye+Sex, ~Hair*Eye)

v1<- lapply(mod.list,
    function(m){cat("\nMODEL FORMULA GIVEN:", paste(m),"\n\n")
                mo<-hllm(m,gmData=gm.hec);
                mof<- fit(mo)
                print(mof)
                }
                )

v2<- lapply(mod.list,
    function(m){cat("\nMODEL FORMULA GIVEN:", paste(m),"\n\n")
                mo<-hllm(m,gmData=gm.hec,marginal=c("Hair","Eye"));
                mof<- fit(mo)
                print(mof)
                mof
                }
                )
