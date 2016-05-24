randLmer <-
function(model,observ) {
  rand_terms<- sapply(findbars(formula(model)),function(x) paste0("(", deparse(x), ")"))
  drop_form<- list()
  rtable<- data.frame(term=rand_terms)
  rtable$d.AIC<- NA;rtable$d.BIC<- NA;rtable$Chi.sq<- NA; rtable$p.value<- NA
if (length(fixef(model)) == 1) {
  resp_term<- sapply(nobars(formula(model)),function(x) deparse(x))[2]
 for (i in 1:length(rand_terms)) {
  drop_form[[i]] <- reformulate(rand_terms[-i],response=resp_term)
if (grepl("\\bREML\\b", summary(model)$methTitle) == T) {
  m_new<- lmer(formula=drop_form[[i]],data=observ)
  p_mod<- anova(model,m_new,refit=F) }
if (grepl("\\bREML\\b", summary(model)$methTitle) == F) {
  m_new<- lmer(formula=drop_form[[i]],data=observ,REML=F)
  p_mod<- anova(model,m_new) }
  rtable[,-1][i,]<-c(p_mod$AIC[1]-p_mod$AIC[2],p_mod$BIC[1]-p_mod$BIC[2],p_mod$Chisq[2],p_mod$'Pr(>Chisq)'[2]) } }
if (length(fixef(model)) > 1) {
  resp_term<- sapply(nobars(formula(model)),function(x) deparse(x))[2]
  fixed_terms<- paste(attributes(terms(model))$term.labels)
 for (i in 1:length(rand_terms)) {
  drop_form[[i]] <- reformulate(c(rand_terms[-i],fixed_terms),response=resp_term)
if (grepl("\\bREML\\b", summary(model)$methTitle) == T) {
  m_new<- lmer(formula=drop_form[[i]],data=observ)
  p_mod<- anova(model,m_new,refit=F) }
if (grepl("\\bREML\\b", summary(model)$methTitle) == F) {
  m_new<- lmer(formula=drop_form[[i]],data=observ,REML=F)
  p_mod<- anova(model,m_new) }
  rtable[,-1][i,]<-c(p_mod$AIC[1]-p_mod$AIC[2],p_mod$BIC[1]-p_mod$BIC[2],p_mod$Chisq[2],p_mod$'Pr(>Chisq)'[2]) } }
  invisible(rtable)
}
