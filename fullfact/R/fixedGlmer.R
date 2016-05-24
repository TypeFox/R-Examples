fixedGlmer <-
function(model,observ,fam_link) {
  rand_terms<- sapply(findbars(formula(model)),function(x) paste0("(", deparse(x), ")"))
  drop_form<- list()
  resp_term<- sapply(nobars(formula(model)),function(x) deparse(x))[2]
  fixed_terms<- paste(attributes(terms(model))$term.labels)
  ftable<- data.frame(term=fixed_terms)
  ftable$d.AIC<- NA;ftable$d.BIC<- NA;ftable$Chi.sq<- NA; ftable$p.value<- NA
 for (i in 1:length(fixed_terms)) {
  drop_form[[i]] <- reformulate(c(rand_terms,fixed_terms[-i]),response=resp_term)
  m_new<- glmer(formula=drop_form[[i]],family=fam_link,data=observ)
  p_mod<- anova(model,m_new)
  ftable[,-1][i,]<-c(p_mod$AIC[1]-p_mod$AIC[2],p_mod$BIC[1]-p_mod$BIC[2],p_mod$Chisq[2],p_mod$'Pr(>Chisq)'[2]) }
  invisible(ftable)
}
