convertBack <-
function(tree,otree,regshifts){
	regs<-repaint(otree,regshifts)
	otree2<-data.frame(as(otree,"data.frame"),regimes=regs,shifts=rep(NA,length(otree@nodes)))
	otree2$shifts[match(names(regshifts),otree2$nodes)]<-1:length(regshifts)
	otree2
}
