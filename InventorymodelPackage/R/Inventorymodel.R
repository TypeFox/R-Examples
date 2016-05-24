Inventorymodel <-
function(model=c("EOQ","EPQ","STI","FOC","MCT","MWHC","MWHC2"),
n=NA,a=NA,av=NA,d=NA,h=NA,m=NA,r=NA,s=NA,K=NA,b=NA,c1=NA,c2=NA,cooperation=c(0,1),allocation=c(0,1)){

if (model=="EOQ"){

print("EOQ")
if (sum(is.na(n)==T|is.na(a)==T|is.na(h)==T)>=1|length(d)!=n|length(h)!=n|(sum(is.na(d))>=1& sum(is.na(m))>=1)){
sol<-c("Error: invalid data")
}else{

if (cooperation==0){
sol<-EOQ(n,a,d,h,m)
      if (allocation==1){
sol<-list(sol,SOC(n,a,d,h,m,model="EOQ",cooperation=0))
names(sol)<-c("*","SOC rule")
}
}
if (cooperation==1){
sol<-EOQcoo(n,a,d,h,m)
if (allocation==1){
sol<-list(sol,SOC(n,a,d,h,m,model="EOQ",cooperation=1))
names(sol)<-c("*","SOC rule")
}
}
}
}
if (model=="EPQ"){

print("EPQ")
if (sum(is.na(n)==T|is.na(a)==T|is.na(h)==T|length(r)!=n|length(s)!=n|length(d)!=n|length(h)!=n|is.na(r)==T|is.na(r)==T)>=1| (sum(is.na(d))>=1& sum(is.na(m))>=1)){
sol<-c("Error: invalid data")
}else{
if (cooperation==0){
sol<-EPQ(n,a,d,h,m,r,s)
if (allocation==1){
sol<-list(sol,SOC(n,a,d,h,m,r,s,model="EPQ",cooperation=0))
names(sol)<-c("*","SOC rule")
}
}
if (cooperation==1){
sol<-EPQcoo(n,a,d,h,m,r,s)
if (allocation==1){
sol<-list(sol,SOC(n,a,d,h,m,r,s,model="EPQ",cooperation=1))
names(sol)<-c("*","SOC rule")
}
}
}}
if (model=="STI"){
print("STI")
if (sum(is.na(n)==T|is.na(a)==T|length(av)!=n|length(d)!=n|length(h)!=n|is.na(h)==T|is.na(av)==T)>=1| (sum(is.na(d))>=1& sum(is.na(m))>=1)){
sol<-c("Error: invalid data")
}else{
if (cooperation==0){sol<-STI(n,a,av,d,h,m)}
if (cooperation==1){sol<-STIcoo(n,a,av,d,h,m)}
if (allocation==1){
print("Optimal solution")
sol<-list(sol,linerulecoalitional(n,a,av,d,h,m))
names(sol)<-c("*","Allocation")
}
}
}
if (model=="FOC"){
print("FOC")
if (sum(is.na(n)==T|is.na(a)==T|is.na(d)==T|is.na(K)==T|length(d)!=n|length(K)!=n)>=1){
sol<-c("Error: invalid data")
}else{
sol<-mfoc(n,a,d,K,cooperation)
if (allocation==1){
sol<-list(sol,shapley_mfoc(n,a,d,K))
names(sol)<-c("","Allocation")
}
}}
if (model=="MCT"){
print("MCT")
if (sum(is.na(n)==T|is.na(a)==T|length(av)!=n|length(d)!=n|length(K)!=n|is.na(av)==T|is.na(d)==T|is.na(K)==T)>=1){
sol<-c("Error: invalid data")
}else{
sol<-mct(n,a,av,d,K,cooperation)
if (allocation==1){
sol<-list(sol,twolines(n,a,av,d,K))
names(sol)<-c("Optimal solution","Allocation two-lines rule")
}
}
}
if (model=="MWHC"){
print("MWHC")
if (sum(is.na(n)==T|is.na(a)==T|is.na(b)==T|length(b)!=n|length(d)!=n|length(K)!=n|is.na(d)==T|is.na(K)==T)>=1){
sol<-c("Error: invalid data")
}else{
sol<-mwhc(n,a,b,d,K,cooperation,allocation)
}
}
if (model=="MWHC2"){
print("MWHC2C")
if (sum(is.na(n)==T|is.na(a)==T|is.na(c1)==T|is.na(c2)==T|is.na(b)==T|length(b)!=n|length(d)!=n|length(K)!=n|is.na(d)==T|is.na(K)==T)>=1){
sol<-c("Error: invalid data")
}else{
sol<-mwhc2c(n,a,b,d,K,c1,c2,cooperation,allocation)
}
}
return(sol)}
