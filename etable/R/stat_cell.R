stat_cell <-
function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min, digits=3, digits2=1){
kwordlist<-c('N','MEAN','SD','MSD', 'VAR','MEDIAN','MD','MAD','MQQ','PROP','RANGE','CV','MAX','MIN','SUM','IQR','MODE','MISS','SKEW','KURT','P1','P2','P2.5','P5','P10','P25','P50','P75','P90','P95','P97.5','P98','P99', 'P99.9','PNM', 'M2SD', 'POP', 'COMB')
if(vars[1]%in%kwordlist) keyword<- vars[1]
if(vars[2]%in%kwordlist) keyword<- vars[2]

if(is.null(x) & !is.null(y)){
x<- y
}



out<-''
if(length(cell_ids)>= n_min){
#########################################
if(keyword=='N'){
if(is.null(w)) w<- rep(1, length(x))
w[which(is.na(x))]<-NA
out<- a.round.ade(sum(w[cell_ids], na.rm=TRUE), 0)
}
#########################################

#########################################
if(keyword=='PROP'){
if(is.null(w)) w<- rep(1, length(x))
w[which(is.na(x))]<-NA
out<- paste(a.round.ade((sum(w[cell_ids], na.rm=TRUE)/sum(w, na.rm=TRUE))*100, digits2), '%', sep='')
}
#########################################

#########################################
if(keyword=='POP'){
if(is.null(w)) w<- rep(1, length(x))
out<- ' '
if(nlevels(as.factor(x))==2){
w[which(is.na(x))]<-NA
x<- as.factor(x)
out<- paste(a.round.ade((sum(w[cell_ids][which(x[cell_ids]==levels(x)[2])], na.rm=TRUE)/(sum(w[cell_ids], na.rm=TRUE)))*100, digits2), '% (', round(sum(w[cell_ids][which(x[cell_ids]==levels(x)[2])], na.rm=TRUE), digits=0), ')', sep='')
}
}
#########################################

#########################################
if(keyword=='COMB'){


if(is.factor(x) & nlevels(x)==2){
if(is.null(w)) w<- rep(1, length(x))
out<- ' '
if(nlevels(as.factor(x))==2){
w[which(is.na(x))]<-NA
x<- as.factor(x)
out<- paste(a.round.ade((sum(w[cell_ids][which(x[cell_ids]==levels(x)[2])], na.rm=TRUE)/(sum(w[cell_ids], na.rm=TRUE)))*100, digits2), '% (', round(sum(w[cell_ids][which(x[cell_ids]==levels(x)[2])], na.rm=TRUE), digits=0), ')', sep='')
}
}


if(is.factor(x) & nlevels(x)>2){
out<- ''
}


if(is.numeric(x)){

if(is.null(w)){
m1<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.5),   na.rm=TRUE, type=8), digits=digits)
q1<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.25),  na.rm=TRUE, type=8), digits=digits)
q3<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.75),  na.rm=TRUE, type=8), digits=digits)
out<- paste(m1,' (',q1,'/',q3,')', sep='')
}

if(!is.null(w)){
m1<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.5)), digits=digits)
q1<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.25)), digits=digits)
q3<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.75)), digits=digits)
out<- paste(m1,' (',q1,'/',q3,')', sep='')
}
}

}
#########################################


#########################################
if(keyword=='RANGE'){
out<- paste(a.format_n.ade(min(as.numeric(x[cell_ids]), na.rm=TRUE),digits),'--',a.format_n.ade(max(as.numeric(x[cell_ids]), na.rm=TRUE),digits))
}
#########################################

#########################################
if(keyword=='MEAN'){
if(is.null(w))   w<- rep(1, length(x))
out<-a.format_n.ade(    wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE),  digits=digits)
}
#########################################

#########################################
if(keyword=='SD'){
if(is.null(w))   w<- rep(1, length(x))
out<- a.format_n.ade(sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE)), digits=digits)
}
#########################################

#########################################
if(keyword=='VAR'){
if(is.null(w))   w<- rep(1, length(x))
out<- a.format_n.ade(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE), digits=digits)
}
#########################################

#########################################
if(keyword=='MEDIAN'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.5),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.5)), digits=digits)
}
}
#########################################

#########################################
if(keyword=='MD'){
if(is.null(w)){
v<- x[cell_ids]
out<-a.format_n.ade(mean(abs(v-mean(v,na.rm=TRUE)), na.rm=TRUE)*1.253, digits=digits)
}
if(!is.null(w)){
v<- x[cell_ids]
vw<-w[cell_ids]
out<-a.format_n.ade(wtd.mean(abs(v-wtd.mean(v,vw)), vw)*1.253, digits=digits)
}
}
#########################################

#########################################
if(keyword=='MAD'){
if(is.null(w)){
v<- x[cell_ids]
out<- a.format_n.ade(mad(v, na.rm = TRUE), digits=digits)
}
if(!is.null(w)){
v<- x[cell_ids]
vw<-w[cell_ids]
out<-a.format_n.ade(wtd.quantile(abs(v-wtd.quantile(v,vw,probs=c(0.5))), vw,probs=c(0.5))*1.4826, digits=digits)
}
}
#########################################

#########################################
if(keyword=='CV'){
if(is.null(w))   w<- rep(1, length(x))
mu  <- wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE)
msd <- sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE))

out<- a.format_n.ade(msd/mu, digits=digits)
}
#########################################

#########################################
if(keyword=='MAX'){
out<- a.format_n.ade(as.numeric(max(x[cell_ids], na.rm=TRUE)),digits)
}
#########################################

#########################################
if(keyword=='MIN'){
out<- a.format_n.ade(as.numeric(min(x[cell_ids], na.rm=TRUE)),digits)
}
#########################################

#########################################
if(keyword=='SUM'){
if(!is.null(w)) out<- a.format_n.ade(wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE)*sum(w[cell_ids][which(!is.na(x[cell_ids]))], na.rm=TRUE),digits)
if(is.null(w))  out<- a.format_n.ade(sum(x[cell_ids], na.rm=TRUE),digits)

}
#########################################


#########################################
#########################################
if(keyword=='P1'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.01),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.01)), digits=digits)
}}
#########################################
if(keyword=='P2'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.02),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.02)), digits=digits)
}}
#########################################
if(keyword=='P2.5'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.025),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.025)), digits=digits)
}}
#########################################
if(keyword=='P5'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.05),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.05)), digits=digits)
}}
#########################################
if(keyword=='P10'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.10),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.10)), digits=digits)
}}
#########################################

#########################################
if(keyword=='P25'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.25),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.25)), digits=digits)
}}
#########################################
#########################################
if(keyword=='P50'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.5),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.5)), digits=digits)
}}
#########################################
if(keyword=='P75'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.75),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.75)), digits=digits)
}}
#########################################
if(keyword=='P90'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.90),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.90)), digits=digits)
}}
#########################################
if(keyword=='P95'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.95),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.95)), digits=digits)
}}
#########################################
if(keyword=='P97.5'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.975),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.975)), digits=digits)
}}

#########################################
if(keyword=='P98'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.98),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.98)), digits=digits)
}}

#########################################
if(keyword=='P99'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.99),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.99)), digits=digits)
}}
#########################################

#########################################
if(keyword=='P99.9'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.999),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.999)), digits=digits)
}}
#########################################

#########################################
if(keyword=='IQR'){
if(is.null(w))  out<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.75),  na.rm=TRUE, type=8)-quantile(x[cell_ids], probs=c(0.25),  na.rm=TRUE, type=8), digits=digits)
if(!is.null(w)){
out<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.75))-wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.25)), digits=digits)
}
}
#########################################

#########################################
if(keyword=='MISS'){
out<- sum(is.na(x[cell_ids]))
}
#########################################

#########################################
if(keyword=='PNM'){
nmiss<- sum(is.na(x[cell_ids]))
nnonm<- sum(!is.na(x[cell_ids]))

out<- paste(a.round.ade((nnonm/(nmiss+nnonm))*100, 1), '%', sep='')

}
#########################################


#########################################
if(keyword=='SKEW'){
if(is.null(w))   w<- rep(1, length(x))
mu  <- wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE)
vsd <- sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE))
n<- sum(!is.na(x[cell_ids]))
skew<- (1/n)*(sum(w[cell_ids]^(3/2)*((x[cell_ids]-mu)/vsd)^3, na.rm=T))
out<-a.format_n.ade(skew, digits=digits)
}
#########################################


#########################################
if(keyword=='KURT'){
if(is.null(w))   w<- rep(1, length(x))
mu  <- wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE)
vsd <- sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE))
n<- sum(!is.na(x[cell_ids]))
kurt<- ((1/n)*(sum((w[cell_ids]^2) * ((x[cell_ids]-mu)/vsd)^4, na.rm=T)))-3
out<-a.format_n.ade(kurt, digits=digits)
}
#########################################

#########################################
if(keyword=='MODE'){
if(is.null(w))      w<- rep(1, length(x))
if(is.factor(x))    v<- factor(x[cell_ids])
if(is.numeric(x))   v<- x[cell_ids]
if(is.character(x)) v<- x[cell_ids]
ul <- sort(unique(v))
tab <- wtd.table(v, w[cell_ids], type='table')
on <- which(tab==max(tab, na.rm=TRUE))
value<- ul[on]
N<- tab[on]
if(is.numeric(x))  out<- paste(a.format_n.ade(value,digits), ' (N:', round(N) , ')',  sep='', collapse='| ')
if(!is.numeric(x)) out<- paste(value, ' (N:', round(N) , ')',  sep='', collapse='|')
}
#########################################


#########################################
if(keyword=='MQQ'){
if(is.null(w)){
m1<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.5),   na.rm=TRUE, type=8), digits=digits)
q1<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.25),  na.rm=TRUE, type=8), digits=digits)
q3<-a.format_n.ade(quantile(x[cell_ids], probs=c(0.75),  na.rm=TRUE, type=8), digits=digits)
out<- paste(m1,' (',q1,'/',q3,')', sep='')
}

if(!is.null(w)){
m1<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.5)), digits=digits)
q1<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.25)), digits=digits)
q3<-a.format_n.ade(wtd.quantile(x[cell_ids], w[cell_ids], probs=c(0.75)), digits=digits)
out<- paste(m1,' (',q1,'/',q3,')', sep='')
}
}
#########################################


#########################################
if(keyword=='MSD'){
if(is.null(w)) w<- rep(1, length(x))
m1<-  a.format_n.ade(    wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE),  digits=digits)
vsd<- a.format_n.ade(sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE)), digits=digits)
out<- paste(m1,' (',vsd,')', sep='')
}
#########################################


#########################################
if(keyword=='M2SD'){
if(is.null(w)) w<- rep(1, length(x))
onesd<- sqrt(wtd.var(x[cell_ids], w[cell_ids], na.rm=TRUE))
amean<- wtd.mean(x[cell_ids], w[cell_ids], na.rm=TRUE)
low<- a.format_n.ade(amean-onesd*2, digits=digits)
upp<- a.format_n.ade(amean+onesd*2, digits=digits)
out<- paste(low,' - ',upp, sep='')
}
#########################################


}

return(out)
}
