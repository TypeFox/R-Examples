require(lmerTest)

# system.time(for(i in 5:19){
#   m <- lmer(TVbo[,i] ~ TVset*Picture+
#             (1|Assessor)+(1|Assessor:TVset) + (1|Assessor:Picture)
#           + (1|Assessor:TVset:Picture), data=TVbo)
#   print(names(TVbo)[i])
#   print(anova(m))  
# 
# })


form <- "  TVset*Picture + (1|Assessor)+(1|Assessor:TVset) + (1|Assessor:Picture) + (1|Assessor:TVset:Picture)"
responses <- names(TVbo)[5:19]
get_models <- function(responses) {
  lapply(responses, function(v) {
    form2 <- as.formula(paste(v, form, sep=" ~ "))    
    lmer(form2, data=TVbo)
  })
}


mlist <- get_models(responses)
if(require(plyr)){
  system.time(ls1 <- llply(mlist, anova, .progress = "text"))
  names(ls1) <- responses
  print(ls1)  
}



