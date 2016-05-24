tabular.ade <-
function(x_vars, xname=NULL, y_vars=NULL, yname=NULL, z_vars=NULL, zname=NULL, rows=NULL, rnames=NULL, cols=NULL, cnames=NULL, w=NULL, data=NULL, FUN, allnames=FALSE, nonames=TRUE, alllabel='Total', inset='?', remove='', n_min=0,...){
if(is.null(xname))    xname<-x_vars
if(is.null(yname))    yname<-y_vars
if(is.null(zname))    zname<-z_vars
if(is.null(rnames))   rnames<-rows[which(rows!='ALL')]
if(is.null(cnames))   cnames<-cols[which(cols!='ALL')]
rowstart<- NULL
colstart<- NULL



########################################
# Fall 1
if(length(x_vars)>1 & length(y_vars)<2){
Big_M<- NULL
rows_sub <-rows[which(rows!='ALL')]
cols_sub <-cols[which(cols!='ALL')]


for(i in 1:length(x_vars)){

sub_M <- tabular.ade(x_vars=x_vars[i], xname=xname[i], y_vars=y_vars, yname=yname, z_vars=z_vars, zname=zname, rows=rows, rnames=rnames, cols=cols, cnames=cnames, w=w, data=data, FUN=FUN, allnames=allnames, nonames=FALSE, alllabel=alllabel, n_min=n_min,...)
if(i==1){
rowstart <- attr(sub_M, 'start cell')[1]
colstart <- attr(sub_M, 'start cell')[2]+1
}
if(i>1 & length(cols_sub)>0)  sub_M <- sub_M[-(1:length(cols_sub)), ]
if(i==1 & !allnames)  sub_M <- cbind(c(rep('', length(cols_sub)), xname[i], rep('', (dim(sub_M)[1]-1-length(cols_sub)))), sub_M)
if(i==1 &  allnames)  sub_M <- cbind(c(rep('', length(cols_sub)), rep(xname[i], (dim(sub_M)[1]-length(cols_sub)))), sub_M)
if(!is.matrix(sub_M)) sub_M<- t(as.matrix(sub_M))
if(i>1 & !allnames)   sub_M <- cbind(c(xname[i], rep('', (dim(sub_M)[1]-1))), sub_M)
if(i>1 &  allnames)   sub_M <- cbind(c(rep(xname[i], (dim(sub_M)[1]))), sub_M)
Big_M <- rbind(Big_M, sub_M)
}

if(nonames){
if(!all(rownames(Big_M)=='')) Big_M<- cbind(rownames(Big_M) ,Big_M)
if(!all(colnames(Big_M)=='')) Big_M<- rbind(colnames(Big_M) ,Big_M)
if(!all(colnames(Big_M)=='')) rowstart<- rowstart+1
if(!all(rownames(Big_M)=='')) colstart<- colstart+1
rownames(Big_M)<- 1:dim(Big_M)[1]
colnames(Big_M)<- 1:dim(Big_M)[2]
}
attributes(Big_M)<- list("dim"=dim(Big_M), 'dimnames'=dimnames(Big_M), 'start cell'=c(rowstart, colstart))
return(Big_M)
}
########################################


########################################
# Fall 2
if(length(y_vars)>1 & length(x_vars)<2){
Big_M<- NULL
rows_sub <-rows[which(rows!='ALL')]
cols_sub <-cols[which(cols!='ALL')]


for(i in 1:length(y_vars)){
sub_M <- tabular.ade(x_vars=x_vars, xname=xname, y_vars=y_vars[i], yname=yname[i], z_vars=z_vars, zname=zname, rows=rows, rnames=rnames, cols=cols, cnames=cnames, w=w, data=data, FUN=FUN, allnames=allnames, nonames=FALSE, alllabel=alllabel, n_min=n_min,...)
if(i==1){
rowstart <- attr(sub_M, 'start cell')[1]+1
colstart <- attr(sub_M, 'start cell')[2]
}
if(i>1 & length(rows_sub)>0)  sub_M <- sub_M[, -(1:length(rows_sub))]
if(i==1 & !allnames) sub_M <- rbind(c(rep('', length(rows_sub)), yname[i], rep('', (dim(sub_M)[2]-1-length(rows_sub)))), sub_M)
if(i==1 & allnames)  sub_M <- rbind(c(rep('', length(rows_sub)), rep(yname[i],, (dim(sub_M)[2]-length(rows_sub)))), sub_M)
if(!is.matrix(sub_M))sub_M<- as.matrix(sub_M)
if(i>1 & !allnames)  sub_M <- rbind(c(yname[i], rep('', (dim(sub_M)[2]-1))), sub_M)
if(i>1 & allnames)   sub_M <- rbind(c(rep(yname[i], (dim(sub_M)[2]))), sub_M)

Big_M <- cbind(Big_M, sub_M)
}

if(nonames){
if(!all(rownames(Big_M)=='')) Big_M<- cbind(rownames(Big_M) ,Big_M)
if(!all(colnames(Big_M)=='')) Big_M<- rbind(colnames(Big_M) ,Big_M)
if(!all(colnames(Big_M)=='')) rowstart<- rowstart+1
if(!all(rownames(Big_M)=='')) colstart<- colstart+1
rownames(Big_M)<- 1:dim(Big_M)[1]
colnames(Big_M)<- 1:dim(Big_M)[2]
}

attributes(Big_M)<- list("dim"=dim(Big_M), 'dimnames'=dimnames(Big_M), 'start cell'=c(rowstart, colstart))
return(Big_M)
}
########################################


########################################
# Fall 3
if(length(x_vars)>1 & length(y_vars)>1){
Big_M   <- NULL
Middle_M<- NULL
rows_sub <-rows[which(rows!='ALL')]
cols_sub <-cols[which(cols!='ALL')]
for(i in 1:length(x_vars)){
for(j in 1:length(y_vars)){

sub_M <- tabular.ade(x_vars=x_vars[i], xname=xname[i], y_vars=y_vars[j], yname=yname[j], z_vars=z_vars, zname=zname, rows=rows, rnames=rnames, cols=cols, cnames=cnames, w=w, data=data, FUN=FUN, allnames=allnames, nonames=FALSE, alllabel=alllabel, n_min=n_min,...)
if(i==1 & j==1){
rowstart <- attr(sub_M, 'start cell')[1]+1
colstart <- attr(sub_M, 'start cell')[2]+1
}

if(j>1 & length(rows_sub)>0)   sub_M <- sub_M[, -(1:length(rows_sub))]
#if(j>1 & length(rows_sub)==0)  sub_M <- sub_M[, -1]
if(j==1 & !allnames) sub_M <- rbind(c(rep('', length(rows_sub)),      yname[j], rep('', (dim(sub_M)[2]-1-length(rows_sub)))), sub_M)
if(j==1 & allnames)  sub_M <- rbind(c(rep('', length(rows_sub)),  rep(yname[j], (dim(sub_M)[2]-length(rows_sub)))), sub_M)
if(!is.matrix(sub_M)) sub_M<- as.matrix(sub_M)
if(j>1 & !allnames)  sub_M <- rbind(c(yname[j], rep('', (dim(sub_M)[2]-1))), sub_M)
if(j>1 & allnames)   sub_M <- rbind(c(rep(yname[j], (dim(sub_M)[2]))), sub_M)

Middle_M <- cbind(Middle_M, sub_M)

}

if(i>1 & length(cols_sub)>0)   Middle_M <- Middle_M[-(1:(length(cols_sub)+1)), ]
if(i>1 & length(cols_sub)==0)  Middle_M <- Middle_M[-1, ]
if(i==1 & !allnames) Middle_M <- cbind(c(rep('', length(cols_sub)+1), xname[i], rep('', (dim(Middle_M)[1]-2-length(cols_sub)))), Middle_M)
if(i==1 & allnames)  Middle_M <- cbind(c(rep('', length(cols_sub)+1), rep(xname[i], (dim(Middle_M)[1]-1-length(cols_sub)))), Middle_M)
if(!is.matrix(Middle_M)) Middle_M<- t(as.matrix(Middle_M))
if(i>1 & !allnames)  Middle_M <- cbind(c(xname[i], rep('', (dim(Middle_M)[1]-1))), Middle_M)
if(i>1 &  allnames)  Middle_M <- cbind(c(rep(xname[i], (dim(Middle_M)[1]))), Middle_M)


Big_M <- rbind(Big_M, Middle_M)
Middle_M <- NULL

}


if(nonames){
if(!all(rownames(Big_M)=='')) Big_M<- cbind(rownames(Big_M) ,Big_M)
if(!all(colnames(Big_M)=='')) Big_M<- rbind(colnames(Big_M) ,Big_M)
if(!all(colnames(Big_M)=='')) rowstart<- rowstart+1
if(!all(rownames(Big_M)=='')) colstart<- colstart+1
rownames(Big_M)<- 1:dim(Big_M)[1]
colnames(Big_M)<- 1:dim(Big_M)[2]
}


attributes(Big_M)<- list("dim"=dim(Big_M), 'dimnames'=dimnames(Big_M), 'start cell'=c(rowstart, colstart))
return(Big_M)
}
########################################




################################################################################
################################################################################
if(length(x_vars)==1 & length(y_vars)<2 & length(z_vars)<2){
vars   <- c(x_vars, y_vars, z_vars)
vnames <- c(xname, yname, zname)
if(is.null(vnames))   vnames<-vars
rows_all<-rows
cols_all<-cols
rows<-rows[which(rows!='ALL')]
cols<-cols[which(cols!='ALL')]

gtr<- function(x, id){
if(id<=nlevels(x))      out<- (x==levels(x)[id])
if(id==(nlevels(x)+1))  out<- (!is.na(x))
return(out)
}


################################################################################
# Dummy erzeugen
nd<- as.factor(rep(1, dim(data)[1]))
rowF1<- rowF2 <-rowF3 <-rowF4 <-rowF5 <-rowF6 <-nd
colF1<- colF2 <-colF3 <-colF4 <-colF5 <-colF6 <-nd
x_all <-y_all <-z_all <-w_all<- NULL
################################################################################


################################################################################
# Variablen abrufen
if(!is.null(rows)){
if(length(rows)>0) eval(parse(text=paste('rowF1<- data$', rows[1], sep='')))
if(length(rows)>1) eval(parse(text=paste('rowF2<- data$', rows[2], sep='')))
if(length(rows)>2) eval(parse(text=paste('rowF3<- data$', rows[3], sep='')))
if(length(rows)>3) eval(parse(text=paste('rowF4<- data$', rows[4], sep='')))
if(length(rows)>4) eval(parse(text=paste('rowF5<- data$', rows[5], sep='')))
if(length(rows)>5) eval(parse(text=paste('rowF6<- data$', rows[6], sep='')))
if(length(rows)>0) rowF1<- as.factor(rowF1)
if(length(rows)>1) rowF2<- as.factor(rowF2)
if(length(rows)>2) rowF3<- as.factor(rowF3)
if(length(rows)>3) rowF4<- as.factor(rowF4)
if(length(rows)>4) rowF5<- as.factor(rowF5)
if(length(rows)>5) rowF6<- as.factor(rowF6)
}
###################
if(!is.null(cols)){
if(length(cols)>0) eval(parse(text=paste('colF1<- data$', cols[1], sep='')))
if(length(cols)>1) eval(parse(text=paste('colF2<- data$', cols[2], sep='')))
if(length(cols)>2) eval(parse(text=paste('colF3<- data$', cols[3], sep='')))
if(length(cols)>3) eval(parse(text=paste('colF4<- data$', cols[4], sep='')))
if(length(cols)>4) eval(parse(text=paste('colF5<- data$', cols[5], sep='')))
if(length(cols)>5) eval(parse(text=paste('colF6<- data$', cols[6], sep='')))

if(length(cols)>0) colF1<- as.factor(colF1)
if(length(cols)>1) colF2<- as.factor(colF2)
if(length(cols)>2) colF3<- as.factor(colF3)
if(length(cols)>3) colF4<- as.factor(colF4)
if(length(cols)>4) colF5<- as.factor(colF5)
if(length(cols)>5) colF6<- as.factor(colF6)
}
###################

Nrows1<-Nrows2<-Nrows3<-Nrows4<-Nrows5<-Nrows6<-1
Ncols1<-Ncols2<-Ncols3<-Ncols4<-Ncols5<-Ncols6<-1

if(nlevels(rowF1)>1){
at<-which(rows_all==rows[1])
Nrows1<- nlevels(rowF1)
if(length(rows_all)>at) if(rows_all[at+1]=='ALL') Nrows1<- Nrows1+1
}
if(nlevels(rowF2)>1){
at<-which(rows_all==rows[2])
Nrows2<- nlevels(rowF2)
if(length(rows_all)>at) if(rows_all[at+1]=='ALL') Nrows2<- Nrows2+1
}
if(nlevels(rowF3)>1){
at<-which(rows_all==rows[3])
Nrows3<- nlevels(rowF3)
if(length(rows_all)>at) if(rows_all[at+1]=='ALL') Nrows3<- Nrows3+1
}
if(nlevels(rowF4)>1){
at<-which(rows_all==rows[4])
Nrows4<- nlevels(rowF4)
if(length(rows_all)>at) if(rows_all[at+1]=='ALL') Nrows4<- Nrows4+1
}
if(nlevels(rowF5)>1){
at<-which(rows_all==rows[5])
Nrows5<- nlevels(rowF5)
if(length(rows_all)>at) if(rows_all[at+1]=='ALL') Nrows5<- Nrows5+1
}
if(nlevels(rowF6)>1){
at<-which(rows_all==rows[6])
Nrows6<- nlevels(rowF6)
if(length(rows_all)>at) if(rows_all[at+1]=='ALL') Nrows6<- Nrows6+1
}


if(nlevels(colF1)>1){
at<-which(cols_all==cols[1])
Ncols1<- nlevels(colF1)
if(length(cols_all)>at) if(cols_all[at+1]=='ALL') Ncols1<- Ncols1+1
}
if(nlevels(colF2)>1){
at<-which(cols_all==cols[2])
Ncols2<- nlevels(colF2)
if(length(cols_all)>at) if(cols_all[at+1]=='ALL') Ncols2<- Ncols2+1
}
if(nlevels(colF3)>1){
at<-which(cols_all==cols[3])
Ncols3<- nlevels(colF3)
if(length(cols_all)>at) if(cols_all[at+1]=='ALL') Ncols3<- Ncols3+1
}
if(nlevels(colF4)>1){
at<-which(cols_all==cols[4])
Ncols4<- nlevels(colF4)
if(length(cols_all)>at) if(cols_all[at+1]=='ALL') Ncols4<- Ncols4+1
}
if(nlevels(colF5)>1){
at<-which(cols_all==cols[5])
Ncols5<- nlevels(colF5)
if(length(cols_all)>at) if(cols_all[at+1]=='ALL') Ncols5<- Ncols5+1
}
if(nlevels(colF6)>1){
at<-which(cols_all==cols[6])
Ncols6<- nlevels(colF6)
if(length(cols_all)>at) if(cols_all[at+1]=='ALL') Ncols6<- Ncols6+1
}

M<- matrix('', Nrows1*Nrows2*Nrows3*Nrows4*Nrows5*Nrows6+length(cols), Ncols1*Ncols2*Ncols3*Ncols4*Ncols5*Ncols6+length(rows))



if(length(vars)>0) eval(parse(text=paste('x_all<- data$', vars[1], sep='')))
if(length(vars)>1) eval(parse(text=paste('y_all<- data$', vars[2], sep='')))
if(length(vars)>2) eval(parse(text=paste('z_all<- data$', vars[3], sep='')))
if(!is.null(w) & is.character(w)) eval(parse(text=paste('w_all<- data$', w, sep='')))
#if(is.null(w))  w_all<- rep(1, length(x_all))                                  #
################################################################################


################################################################################
# Namengeben Rows
rrun_1<-length(cols)+1
save_r_1<-save_r_2<-save_r_3<-save_r_4<-save_r_5<-save_r_6<-'no_factor_ad'

# Namengeben cols
crun_1<-length(rows)+1
save_c_1<-save_c_2<-save_c_3<-save_c_4<-save_c_5<-save_c_6<-'no_factor_ad'
################################################################################



################################################################################
# Berechnung                                                                   #
rowid<-length(cols)
colid<-length(rows)
for(i in 1:Nrows1){
for(k in 1:Nrows2){
for(l in 1:Nrows3){
for(i2 in 1:Nrows4){
for(k2 in 1:Nrows5){
for(l2 in 1:Nrows6){

rowid <- rowid+1

#######################
# Rows subnames
if(!allnames){
if(Nrows1>1){
if(i <= nlevels(rowF1)){
if(levels(rowF1)[i]!=save_r_1) M[rrun_1, 1] <- levels(rowF1)[i]
save_r_1     <- levels(rowF1)[i]
}
if(i == (nlevels(rowF1)+1)){
if(alllabel!=save_r_1) M[rrun_1, 1] <- alllabel
save_r_1     <- alllabel
}
}
if(Nrows2>1){
if(k <= nlevels(rowF2)){
if(levels(rowF2)[k]!=save_r_2) M[rrun_1, 2] <- levels(rowF2)[k]
save_r_2     <- levels(rowF2)[k]
}
if(k == (nlevels(rowF2)+1)){
if(alllabel!=save_r_2) M[rrun_1, 2] <- alllabel
save_r_2     <- alllabel
}
}
if(Nrows3>1){
if(l <= nlevels(rowF3)){
if(levels(rowF3)[l]!=save_r_3) M[rrun_1, 3] <- levels(rowF3)[l]
save_r_3     <- levels(rowF3)[l]
}
if(l == (nlevels(rowF3)+1)){
if(alllabel!=save_r_3) M[rrun_1, 3] <- alllabel
save_r_3     <- alllabel
}
}
if(Nrows4>1){
if(i2 <= nlevels(rowF4)){
if(levels(rowF4)[i2]!=save_r_4) M[rrun_1, 4] <- levels(rowF4)[i2]
save_r_4     <- levels(rowF4)[i2]
}
if(i2 == (nlevels(rowF4)+1)){
if(alllabel!=save_r_4) M[rrun_1, 4] <- alllabel
save_r_4     <- alllabel
}
}
if(Nrows5>1){
if(k2 <= nlevels(rowF5)){
if(levels(rowF5)[k2]!=save_r_5) M[rrun_1, 5] <- levels(rowF5)[k2]
save_r_5     <- levels(rowF5)[k2]
}
if(k2 == (nlevels(rowF5)+1)){
if(alllabel!=save_r_5) M[rrun_1, 5] <- alllabel
save_r_5     <- alllabel
}
}

if(Nrows6>1){
if(l2 <= nlevels(rowF6)){
if(levels(rowF6)[l2]!=save_r_6) M[rrun_1, 6] <- levels(rowF6)[l2]
save_r_6     <- levels(rowF6)[l2]
}
if(l2 == (nlevels(rowF6)+1)){
if(alllabel!=save_r_6) M[rrun_1, 6] <- alllabel
save_r_6     <- alllabel
}
}
}

if(allnames){
if(Nrows1>1){
if(i <= nlevels(rowF1))     M[rrun_1, 1] <- levels(rowF1)[i]
if(i == (nlevels(rowF1)+1)) M[rrun_1, 1] <- alllabel
}
if(Nrows2>1){
if(k <= nlevels(rowF2))     M[rrun_1, 2] <- levels(rowF2)[k]
if(k == (nlevels(rowF2)+1)) M[rrun_1, 2] <- alllabel
}
if(Nrows3>1){
if(l <= nlevels(rowF3))     M[rrun_1, 3] <- levels(rowF3)[l]
if(l == (nlevels(rowF3)+1)) M[rrun_1, 3] <- alllabel
}
if(Nrows4>1){
if(i2 <= nlevels(rowF4))     M[rrun_1, 4] <- levels(rowF4)[i2]
if(i2 == (nlevels(rowF4)+1)) M[rrun_1, 4] <- alllabel
}
if(Nrows5>1){
if(k2 <= nlevels(rowF5))     M[rrun_1, 5] <- levels(rowF5)[k2]
if(k2 == (nlevels(rowF5)+1)) M[rrun_1, 5] <- alllabel
}
if(Nrows6>1){
if(l2 <= nlevels(rowF6))     M[rrun_1, 6] <- levels(rowF6)[l2]
if(l2 == (nlevels(rowF6)+1)) M[rrun_1, 6] <- alllabel
}
}

rrun_1<- rrun_1+1
#######################

#####################################################
for(m in 1:Ncols1){
for(n in 1:Ncols2){
for(o in 1:Ncols3){
for(m2 in 1:Ncols4){
for(n2 in 1:Ncols5){
for(o2 in 1:Ncols6){
colid <- colid+1



#######################
# Cols subnames
if(i==1 & k==1 & l==1 & i2==1 & k2==1 & l2==1 ){
if(!allnames){
if(Ncols1>1){
if(m <= nlevels(colF1)){
if(levels(colF1)[m]!=save_c_1) M[1, crun_1] <- levels(colF1)[m]
save_c_1     <- levels(colF1)[m]
}
if(m == (nlevels(colF1)+1)){
if(alllabel!=save_c_1)  M[1, crun_1] <- alllabel
save_c_1     <- alllabel
}
}
if(Ncols2>1){
if(n <= nlevels(colF2)){
if(levels(colF2)[n]!=save_c_2) M[2, crun_1] <- levels(colF2)[n]
save_c_2     <- levels(colF2)[n]
}
if(n == (nlevels(colF2)+1)){
if(alllabel!=save_c_2)  M[2, crun_1] <- alllabel
save_c_2     <- alllabel
}
}
if(Ncols3>1){
if(o <= nlevels(colF3)){
if(levels(colF3)[o]!=save_c_3) M[3, crun_1] <- levels(colF3)[o]
save_c_3     <- levels(colF3)[o]
}
if(o == (nlevels(colF3)+1)){
if(alllabel!=save_c_3)  M[3, crun_1] <- alllabel
save_c_3     <- alllabel
}
}
if(Ncols4>1){
if(m2 <= nlevels(colF4)){
if(levels(colF4)[m2]!=save_c_4) M[4, crun_1] <- levels(colF4)[m2]
save_c_4     <- levels(colF4)[m2]
}
if(m2 == (nlevels(colF4)+1)){
if(alllabel!=save_c_4)  M[4, crun_1] <- alllabel
save_c_4     <- alllabel
}
}
if(Ncols5>1){
if(n2 <= nlevels(colF5)){
if(levels(colF5)[n2]!=save_c_5) M[5, crun_1] <- levels(colF5)[n2]
save_c_5     <- levels(colF5)[n2]
}
if(n2 == (nlevels(colF5)+1)){
if(alllabel!=save_c_5)  M[5, crun_1] <- alllabel
save_c_5     <- alllabel
}
}
if(Ncols6>1){
if(o2 <= nlevels(colF6)){
if(levels(colF6)[o2]!=save_c_6) M[6, crun_1] <- levels(colF6)[o2]
save_c_6     <- levels(colF6)[o2]
}
if(o2 == (nlevels(colF6)+1)){
if(alllabel!=save_c_6)  M[6, crun_1] <- alllabel
save_c_6     <- alllabel
}
}
}


if(allnames){
if(Ncols1>1){
if(m <= nlevels(colF1))     M[1, crun_1] <- levels(colF1)[m]
if(m == (nlevels(colF1)+1)) M[1, crun_1] <- alllabel
}
if(Ncols2>1){
if(n <= nlevels(colF2))     M[2, crun_1] <- levels(colF2)[n]
if(n == (nlevels(colF2)+1)) M[2, crun_1] <- alllabel
}
if(Ncols3>1){
if(o <= nlevels(colF3))     M[3, crun_1] <- levels(colF3)[o]
if(o == (nlevels(colF3)+1)) M[3, crun_1] <- alllabel
}
if(Ncols4>1){
if(m2 <= nlevels(colF4))     M[4, crun_1] <- levels(colF4)[m2]
if(m2 == (nlevels(colF4)+1)) M[4, crun_1] <- alllabel
}
if(Ncols5>1){
if(n2 <= nlevels(colF5))     M[5, crun_1] <- levels(colF5)[n2]
if(n2 == (nlevels(colF5)+1)) M[5, crun_1] <- alllabel
}
if(Ncols6>1){
if(o2 <= nlevels(colF6))     M[6, crun_1] <- levels(colF6)[o2]
if(o2 == (nlevels(colF6)+1)) M[6, crun_1] <- alllabel
}
}


crun_1<- crun_1+1
}
#######################

#####################################################
#######################
cell_ids  <-which(gtr(rowF1,i) & gtr(rowF2,k) & gtr(rowF3,l) & gtr(rowF4,i2) & gtr(rowF5,k2) & gtr(rowF6,l2) & gtr(colF1,m) & gtr(colF2,n) & gtr(colF3,o) & gtr(colF4,m2) & gtr(colF5,n2) & gtr(colF6,o2))

row_ids   <-which(gtr(rowF1,i) & gtr(rowF2,k) & gtr(rowF3,l) & gtr(rowF4,i2) & gtr(rowF5,k2) & gtr(rowF6,l2))
col_ids   <-which(gtr(colF1,m) & gtr(colF2,n) & gtr(colF3,o) & gtr(colF4,m2) & gtr(colF5,n2) & gtr(colF6,o2))

if(is.null(rowstart)) rowstart<- rowid
if(is.null(colstart)) colstart<- colid

M[rowid, colid]<- FUN(x=x_all, y=y_all, z=z_all, w=w_all, cell_ids=cell_ids, row_ids=row_ids, col_ids=col_ids, vnames=vnames, vars=vars, n_min=n_min, ...)

if(nchar(remove)>0) M[rowid, colid]<- gsub(paste('[', remove, ']', sep=''), '', M[rowid, colid])
if(nchar(inset)>1)  M[rowid, colid]<- gsub('[?]', M[rowid, colid], inset)


#####################################################

}}}}}}
colid<-length(rows)
}}}}}}
################################################################################


rownames(M)<-rep('', dim(M)[1])
colnames(M)<-rep('', dim(M)[2])
if(!is.null(cols)) rownames(M)[1:length(cols)]<- cnames
if(!is.null(rows)) colnames(M)[1:length(rows)]<- rnames

if(nonames){
if(!all(rownames(M)=='')) M<- cbind(rownames(M) ,M)
if(!all(colnames(M)=='')) M<- rbind(colnames(M) ,M)
if(!all(colnames(M)==''))  rowstart<- rowstart+1
if(!all(rownames(M)==''))  colstart<- colstart+1
rownames(M)<- 1:dim(M)[1]
colnames(M)<- 1:dim(M)[2]
}


attributes(M)<- list("dim"=dim(M), 'dimnames'=dimnames(M), 'start cell'=c(rowstart, colstart))
return(M)

}

}
