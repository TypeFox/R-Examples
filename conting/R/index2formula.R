index2formula <-
function(index,maximal.mod){
terms<-attributes(maximal.mod$x)$assign+1					## terms in maximal model

indo.terms<-unique(terms[index==1])					## terms in formula model

modA<-maximal.mod

term.labels<-attr(summary(modA)$terms,"term.labels")[indo.terms[-1]-1]	## term labels in formula model

form<-paste0("y~",term.labels[1])					## cobble together the formula
if(length(term.labels)>1){
for(i in 2:length(term.labels)){
form<-paste0(form,"+",term.labels[i])}}
form<-as.formula(form)

form}
