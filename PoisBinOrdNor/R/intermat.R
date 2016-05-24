intermat<-function(no_pois,no_bin,no_ord,no_norm,corr_mat=NULL,prop_vec_bin=NULL,prop_vec_ord=NULL,lam_vec=NULL,nor_mean=NULL,nor_var=NULL){

if(no_pois !=0) {
n1=no_pois; lambdavec=lam_vec} 
if(no_bin  !=0) {
n2=no_bin; p1=prop_vec_bin}
if(no_ord  !=0) {
n3=no_ord; p2=prop_vec_ord}
if(no_norm !=0) {
n4=no_norm; normean=nor_mean; norvar=nor_var}

d=n1+n2+n3+n4

if(n1==0 && n2==0 &&n3==0 &&n4==0){stop("Number of variables cannot all be zero!\n")}

validation_specs(n1, n2, n3, n4, corr_mat, p1, p2, lambdavec, normean, norvar)


inter.mat=diag(nrow(corr_mat))

if(n1!=0 && n2!=0 && n3!=0 && n4!=0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}

#bo
 inter.mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)]=ordcont(c(as.list(p1), as.list(p2)), corr_mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)])$SigmaC

#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]

#pb
 for (i in (n1+1):(n1+n2)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p1[i-n1], corr_mat[i,j])}}}

#po
 for (i in (n1+n2+1):(n1+n2+n3)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p2[[i-n1-n2]], corr_mat[i,j])}}}

#pn
 for (i in (n1+n2+n3+1):d){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pn(lambdavec[j], corr_mat[i,j])}}} 

#bn
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+1):(n1+n2)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bn(p1[j-n1], corr_mat[i,j])}}}

#on
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+n2+1):(n1+n2+n3)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4on(p2[[j-n1-n2]], corr_mat[i,j])}}}
}

if(n1!=0 && n2!=0 && n3!=0 && n4==0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}

#bo
 inter.mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)]=ordcont(c(as.list(p1), as.list(p2)), corr_mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)])$SigmaC

#pb
 for (i in (n1+1):(n1+n2)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p1[i-n1], corr_mat[i,j])}}}

#po
 for (i in (n1+n2+1):(n1+n2+n3)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p2[[i-n1-n2]], corr_mat[i,j])}}}
}

if(n1!=0 && n2!=0 && n3==0 && n4!=0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}

#bb
 for (i in (n1+1):(n1+n2)){
 for (j in (n1+1):(n1+n2)){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bb(p1[i-n1], p1[j-n1], corr_mat[i,j])}}}

#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]

#pb
 for (i in (n1+1):(n1+n2)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p1[i-n1], corr_mat[i,j])}}}

#pn
 for (i in (n1+n2+n3+1):d){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pn(lambdavec[j], corr_mat[i,j])}}} 

#bn
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+1):(n1+n2)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bn(p1[j-n1], corr_mat[i,j])}}}
}

if(n1!=0 && n2!=0 && n3==0 && n4==0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}

#bb
 for (i in (n1+1):(n1+n2)){
 for (j in (n1+1):(n1+n2)){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bb(p1[i-n1], p1[j-n1], corr_mat[i,j])}}}

#pb
 for (i in (n1+1):(n1+n2)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p1[i-n1], corr_mat[i,j])}}}
}


if(n1!=0 && n2==0 && n3!=0 && n4!=0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}

#oo
inter.mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)]=ordcont(as.list(p2), corr_mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)])$SigmaC

#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]

#po
 for (i in (n1+n2+1):(n1+n2+n3)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p2[[i-n1-n2]], corr_mat[i,j])}}}

#pn
 for (i in (n1+n2+n3+1):d){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pn(lambdavec[j], corr_mat[i,j])}}} 

#on
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+n2+1):(n1+n2+n3)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4on(p2[[j-n1-n2]], corr_mat[i,j])}}}
}

if(n1!=0 && n2==0 && n3!=0 && n4==0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}

#oo
inter.mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)]=ordcont(as.list(p2), corr_mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)])$SigmaC

#po
 for (i in (n1+n2+1):(n1+n2+n3)){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pbo(lambdavec[j], p2[[i-n1-n2]], corr_mat[i,j])}}}
}

if(n1!=0 && n2==0 && n3==0 && n4!=0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}

#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]

#pn
 for (i in (n1+n2+n3+1):d){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pn(lambdavec[j], corr_mat[i,j])}}} 
}

if(n1!=0 && n2==0 && n3==0 && n4==0) {
#pp
 for (i in 1:n1){
 for (j in 1:n1){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4pp(lambdavec[i], lambdavec[j], corr_mat[i,j])}}}
}


if(n1==0 && n2!=0 && n3!=0 && n4!=0) {
#bo
 inter.mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)]=ordcont(c(as.list(p1), as.list(p2)), corr_mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)])$SigmaC

#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]

#bn
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+1):(n1+n2)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bn(p1[j-n1], corr_mat[i,j])}}}

#on
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+n2+1):(n1+n2+n3)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4on(p2[[j-n1-n2]], corr_mat[i,j])}}}
}

if(n1==0 && n2!=0 && n3!=0 && n4==0) {
#bo
 inter.mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)]=ordcont(c(as.list(p1), as.list(p2)), corr_mat[(n1+1):(n1+n2+n3), (n1+1):(n1+n2+n3)])$SigmaC
}

if(n1==0 && n2!=0 && n3==0 && n4!=0) {
#bb
 for (i in (n1+1):(n1+n2)){
 for (j in (n1+1):(n1+n2)){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bb(p1[i-n1], p1[j-n1], corr_mat[i,j])}}}

#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]

#bn
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+1):(n1+n2)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bn(p1[j-n1], corr_mat[i,j])}}}
}

if(n1==0 && n2!=0 && n3==0 && n4==0) {
#bb
 for (i in (n1+1):(n1+n2)){
 for (j in (n1+1):(n1+n2)){          
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4bb(p1[i-n1], p1[j-n1], corr_mat[i,j])}}}
}

if(n1==0 && n2==0 && n3!=0 && n4!=0) {
#oo
inter.mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)]=ordcont(as.list(p2), corr_mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)])$SigmaC

#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]

#on
 for (i in (n1+n2+n3+1):d){          
 for (j in (n1+n2+1):(n1+n2+n3)){
 if (i!=j){inter.mat[i,j]=inter.mat[j,i]=corr.nn4on(p2[[j-n1-n2]], corr_mat[i,j])}}}
}

if(n1==0 && n2==0 && n3!=0 && n4==0) {
#oo
inter.mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)]=ordcont(as.list(p2), corr_mat[(n1+n2+1):(n1+n2+n3), (n1+n2+1):(n1+n2+n3)])$SigmaC
}

if(n1==0 && n2==0 && n3==0 && n4!=0) {
#nn
 inter.mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]=corr_mat[(n1+n2+n3+1):d, (n1+n2+n3+1):d]
}


if(!is.positive.definite(inter.mat)){
warning( "Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
inter.mat=as.matrix(nearPD(inter.mat, corr = TRUE, keepDiag = TRUE)$mat)
inter.mat = (inter.mat+t(inter.mat))/2
}
return(inter.mat)
}

