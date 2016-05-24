SpaceFilling <-
function(asch){

fun1<-function() {
n<-readline("Number of lines of association schemes array :\n")
l<-readline("Number of columns of association schemes array :\n")
n<-as.integer(n);l<-as.integer(l)
return(c(n,l))}

fun2<-function() {
n<-readline("Number of lines of association schemes array :\n")
l<-readline("Number of columns of association schemes array :\n")
w<-readline("Number of the association scheme arrays :\n")
n<-as.integer(n);l<-as.integer(l);w<-as.integer(w)
return(c(n,l,w))}

# Si Div
if (asch == "Div"){
V<-fun1();n<-V[1];l<-V[2]
s<-n*l;A<-matrix(1:s, ncol = V[2], byrow=TRUE)
SF<-matrix(ncol=s,nrow=s)
for (d in 1:s) {
SF[d,d]<-1
for (dd in 1:s){
D<-which(A==d); d1<-D%%n ; if (d1==0){d1<-n};DD<-which(A==dd); d2<-DD%%n ; if (d2==0){d2<-n}
if (d1==d2) {SF[d,dd]<-1;SF[dd,d]<-1}
else {SF[d,dd]<-2;SF[dd,d]<-2}}}}

##SI Rect
if (asch == "Rect"){
V<-fun1();n<-V[1];l<-V[2];s<-n*l;A<-matrix(1:s, ncol =l, byrow=TRUE)
SF<-matrix(ncol=s,nrow=s)
for (d in 1:s) {
SF[d,d]<-1
for (dd in 1:s){
B<-t(A)
D<-which(A==d); d1<-D%%n ; if (d1==0){d1<-n};DD<-which(A==dd); d2<-DD%%n 
if (d2==0){d2<-n}
D1<-which(B==d); d11<-D1%%l ; if (d11==0){d11<-l};DD1<-which(B==dd); d12<-DD1%%l
if (d12==0){d12<-l}
if (d1==d2) {SF[d,dd]<-1;SF[dd,d]<-1}
if (is.na(SF[d,dd])==TRUE){
if (d11==d12) {SF[d,dd]<-2;SF[dd,d]<-2}}}}
for (d in 1:s) {
for (dd in 1:s){
if (is.na(SF[d,dd])==TRUE){
SF[d,dd]<-3;SF[dd,d]<-3}}}}

##si Nestdiv 
if (asch == "Nestdiv"){
V<-fun2();n<-V[1];l<-V[2];w<-V[3]
s<-l*n;A<-NULL;S<-l*n*w
SF<-matrix(ncol=S,nrow=S)
for (i in 1:w){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE)
z<-(i-1)*s
A[[i]]<-A[[i]]+z};B<-Reduce("rbind",A)
for (i in 1:w) {
a<-A[[i]];mi<-min(a);ma<-max(a)
for (d in mi:ma){
for (dd in mi:ma){
D<-which(a==d); d1<-D%%n ; if (d1==0){d1<-n}
DD<-which(a==dd); d2<-DD%%n ; if (d2==0){d2<-n}
if (d1==d2) {SF[d,dd]<-1;SF[dd,d]<-1}
else {SF[d,dd]<-2;SF[dd,d]<-2}}}}
for (d in 1:S) {
for (dd in 1:S){
if (is.na(SF[d,dd])==TRUE){
SF[d,dd]<-3;SF[dd,d]<-3}}}}


#### SI RightAng
if (asch == "RightAng"){
V<-fun2();n<-V[1];l<-V[2];w<-V[3];s<-l*n;A<-NULL;S<-l*n*w
SF<-matrix(ncol=S,nrow=S)
for (i in 1:w){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE)
z<-(i-1)*s
A[[i]]<-A[[i]]+z};B<-Reduce("rbind",A)
for (i in 1:w) {
a<-A[[i]];mi<-min(a);ma<-max(a)
for (d in mi:ma){
for (dd in mi:ma){
D<-which(a==d); d1<-D%%n ; if (d1==0){d1<-n}
DD<-which(a==dd); d2<-DD%%n ; if (d2==0){d2<-n}
if (d1==d2) {SF[d,dd]<-1;SF[dd,d]<-1}
else {SF[d,dd]<-2;SF[dd,d]<-2}}
for (i in 1:w) {
if (i < w){
b<-A[[i+1]];mib<-min(b);mab<-max(b)
for (db in mib:mab){
DB<-which(b==db); db2<-DB%%n ; if (db2==0){db2<-n}
if (d1==db2) {if (is.na(SF[d,db])==TRUE){
SF[d,db]<-3;SF[db,d]<-3}}
else {if (is.na(SF[d,db])==TRUE){
SF[d,db]<-4;SF[db,d]<-4}}}}}}}}


#### SI GrectRightAng4 
if (asch == "GrectRightAng4"){
V<-fun2();n<-V[1];l<-V[2];w<-V[3];s<-l*n;A<-NULL;S<-l*n*w;SF<-matrix(ncol=S,nrow=S)
for (i in 1:w){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE);z<-(i-1)*s
A[[i]]<-A[[i]]+z};B<-Reduce("rbind",A)
for (i in 1:w) {
a<-A[[i]];mi<-min(a);ma<-max(a)
B<-t(a)
for (d in mi:ma){
for (dd in mi:ma){
D<-which(a==d); d1<-D%%n ; if (d1==0){d1<-n}
DD<-which(a==dd); d2<-DD%%n ; if (d2==0){d2<-n}
if (d1==d2) {SF[d,dd]<-1;SF[dd,d]<-1}
D1<-which(B==d); d11<-D1%%n ; if (d11==0){d11<-n}
DD1<-which(B==dd); d21<-DD1%%n ; if (d21==0){d21<-n}
if (d11==d21) 
if (is.na(SF[d,dd])==TRUE){
{SF[d,dd]<-2;SF[dd,d]<-2}}
if (is.na(SF[d,dd])==TRUE){
SF[d,dd]<-3;SF[dd,d]<-3}}}}
for (d in 1:S) {
for (dd in 1:S){
if (is.na(SF[d,dd])==TRUE){
SF[d,dd]<-4;SF[dd,d]<-4}}}}

## SI GrectRightAng5 
if (asch == "GrectRightAng5"){
V<-fun2();n<-V[1];l<-V[2];w<-V[3];s<-l*n;A<-NULL;S<-l*n*w;SF<-matrix(ncol=S,nrow=S)
for (i in 1:w){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE);z<-(i-1)*s
A[[i]]<-A[[i]]+z};B<-Reduce("rbind",A);SF<-matrix(ncol=S,nrow=S)
for (i in 1:w) {
a<-A[[i]];mi<-min(a);ma<-max(a);B<-t(a)
for (d in mi:ma){
for (dd in mi:ma){
D<-which(a==d); d1<-D%%n ; if (d1==0){d1<-n}
DD<-which(a==dd); d2<-DD%%n ; if (d2==0){d2<-n}
if (d1==d2) {SF[d,dd]<-1;SF[dd,d]<-1}
D1<-which(B==d); d11<-D1%%n ; if (d11==0){d11<-n}
DD1<-which(B==dd); d21<-DD1%%n ; if (d21==0){d21<-n}
if (d11==d21) 
if (is.na(SF[d,dd])==TRUE){
{SF[d,dd]<-2;SF[dd,d]<-2}}
if (is.na(SF[d,dd])==TRUE){
SF[d,dd]<-3;SF[dd,d]<-3}}
for (i in 1:w) {
if (i < w){
bb<-A[[i+1]];mib<-min(bb);mab<-max(bb)
for (db in mib:mab){
DB<-which(bb==db); db2<-DB%%n ; if (db2==0){db2<-n}
if (d1==db2) {if (is.na(SF[d,db])==TRUE){SF[d,db]<-4;SF[db,d]<-4}}
else {if (is.na(SF[d,db])==TRUE){SF[d,db]<-5;SF[db,d]<-5}}}}}}}}


## SI GrectRightAng7 
if (asch == "GrectRightAng7"){
V<-fun2();n<-V[1];l<-V[2];w<-V[3];s<-l*n;S<-l*n*w;A<-NULL
for (i in 1:w){
A[[i]]<-matrix(1:s, ncol=l, byrow=TRUE);z<-(i-1)*s
A[[i]]<-A[[i]]+z};B<-Reduce("rbind",A);SF<-matrix(ncol=S,nrow=S)
for (i in 1:w) {
a<-A[[i]];mi<-min(a);ma<-max(a);B<-t(a)
for (d in mi:ma){
for (dd in mi:ma){
D<-which(a==d); d1<-D%%n ; if (d1==0){d1<-n};DD<-which(a==dd); d2<-DD%%n
if (d2==0){d2<-n}
if (d1==d2) {SF[d,dd]<-1;SF[dd,d]<-1}
D1<-which(B==d); d11<-D1%%n ; if (d11==0){d11<-n};DD1<-which(B==dd);d21<-DD1%%n
if (d21==0){d21<-n}
if (d11==d21) 
if (is.na(SF[d,dd])==TRUE){
{SF[d,dd]<-2;SF[dd,d]<-2}}
if (is.na(SF[d,dd])==TRUE){
SF[d,dd]<-3;SF[dd,d]<-3}}
for (i in 1:w) {
if (i < w){
bb<-A[[i+1]];mib<-min(bb);mab<-max(bb)
for (db in mib:mab){
B2<-t(bb);DB<-which(bb==db); db2<-DB%%n ; if (db2==0){db2<-n}
n1<-which(B2==db); n21<-n1%%n ; if (n21==0){n21<-n}
if (D==DB) {if (is.na(SF[d,db])==TRUE){SF[d,db]<-4;SF[db,d]<-4}}
if (d11==db2) {if (is.na(SF[d,db])==TRUE){SF[d,db]<-5;SF[db,d]<-5}}
if (d1==n21) 
if (is.na(SF[d,db])==TRUE){
{SF[d,db]<-6;SF[db,d]<-6}}}}}}}
for (d in 1:S) {
for (dd in 1:S){
if (is.na(SF[d,dd])==TRUE){
SF[d,dd]<-7;SF[dd,d]<-7}}}}
NN<-max(SF)
RR<-dim(SF)[1]
return(list(SFDesign=SF,Runs=RR,Factors=RR,Levels=NN))}

