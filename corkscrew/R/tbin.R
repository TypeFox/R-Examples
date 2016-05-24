tbin <-
function(dv,idv,train,min.obs, plot.bin = c(TRUE,FALSE),method = c(1,2,3)) {                
    nr = function(vec,txt){
if (sum(which(vec == txt))!=0) {
                vec = vec[-which(vec == txt)]
return(vec)
      }else {
    vec
  }
  
}

                if(length(plot.bin)>1){
    plot.bin = FALSE
 }
if(length(method)>1){
    method = 1
 }

 # Determining bins using t-test, when bin2 function gives more than 32 variables
 xyz.bin3 = function(x,method){
     
 tmat = matrix(rep(0,length(levels(x))*length(levels(x))),length(levels(x)),length(levels(x)))
 name = levels(x)
 for (i in 1:length(name)){
 for (j in 1:length(name)) {
 if (i!=j) {
    if (length(train[which(x==name[i]),dv])<5000){
     s1 = shapiro.test(train[which(x==name[i]),dv]) 
 pval1 = s1$p.value
 } else {
     pval1 = 1
 }
 if ( length(train[which(x==name[j]),dv])<5000 ){
 s2 = shapiro.test(train[which(x==name[j]),dv]) 
 pval2 = s2$p.value
 } else {
     pval2 = 1
 }
 
        if (pval1>0.05 & pval2>0.05){
 test = t.test(train[which(x==name[i]),dv],train[which(x==name[j]),dv],alternative = 
"two.sided",paired = FALSE)
 } else {
 test = wilcox.test(train[which(x==name[i]),dv],train[which(x==name[j]),dv],paired = FALSE)
 }
 if (test$p.value<=0.05){
 tmat[i,j] = 0 
 } else {
 tmat[i,j] = 1
 }
 } else {
 tmat[i,j] = 1
 }
 } 
 }

 colnames(tmat) = levels(x)
 rownames(tmat) = levels(x)
 
 # leaving levels small & Other from combining 
 if(sum(name %in% "small")>0){
 tmat[,"small"] = 0
 tmat["small",] = 0
 }
 if(sum(name %in% "others")>0){
 tmat[,"others"] = 0
 tmat["others",] = 0
 }
 g <- graph.adjacency(tmat, mode = "undirected")
 g = simplify(g)
 if(method == 1){
fc <- fastgreedy.community(g)
 } else if (method == 2){
fc <- walktrap.community(g)
 } else if (method == 3) {
fc <- edge.betweenness.community(g)
 }
 com<-cutat(fc, steps= which.max(fc$modularity)-1)
 bucket = membership(fc)
 
 name = nr(name,"small")
 name = nr(name,"others")
 
 for (i in 1:length(bucket)) {

 levels(x)[levels(x)==name[i]] <- paste("cat",bucket[i],sep = "_")
 
 }

return(x) 
}

# removing variables with low sample sizes 
xyz.bin4 = function(x,y){
if (is.factor(x)) {
 t = table(x)
 name = levels(x)
 s = rep(0,length(t))
 for (i in 1:length(t)) {
 if(t[i] <= y){
 s[i] = 1 
 } } 
name = name[s==1]
for (j in 1:length(name)) {
levels(x)[levels(x)==name[j]] <- "small"
 }
return(x) 
} else {
return(x)
}

}

    

     heat = function(x,y,z,train,temp,ptitle) {
       
chick = data.frame(cbind(train[x],train[y],temp[z]))
chick = chick[order(chick[,3]),]
rn = levels(chick[,3])
n = length(rn)

cn = levels(chick[,2])
m = length(cn)

mat = matrix(rep(0,m*n),n,m)

name = names(chick)
measurevar <- as.character(name[1])
groupvars  <- c(as.character(name[3]))
form = as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))

 mval = aggregate(form,chick,mean)
 mval[,2] = round(mval[,2],0)
 rn_new = rn
for(i in 1:length(rn)){
   ind = which(mval[,1]== rn[i])
   rn_new[i] = paste(sprintf("%s(%i)",rn[i],mval[ind,2]))
 }
rownames(mat) = rn_new
colnames(mat) = cn
 
                             
 
 
for(i in 1:n){
 for (j in 1:m) {
 ind_r = which(chick[,3]==rn[i])
 ind_c = which(chick[,2] == cn[j])
 ind = intersect(ind_r,ind_c)
 mat[i,j] = mean(chick[ind,1])
 if (is.nan(mat[i,j]) | is.na(mat[i,j])){
 mat[i,j] = 0
 }
 }
 }
 
 mat1 = t(mat)
     mat = t(mat1[do.call(order,as.data.frame(mat1)),])
 mat = mat[do.call(order,as.data.frame(mat)),]

my_palette <- colorRampPalette(c("white", "yellow", "green"))(n = 500)
heatmap.2(mat,
  main = sprintf("Binning for %s",ptitle), 
  Colv = FALSE,
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins =c(12,12),     
  col=my_palette,       
  dendrogram = c("none"))            
  
}


 # Try Catch for errors arising using t-test binning
     heat.map <- function(dv,x,y,train,temp,ptitle) (
  out <- tryCatch(
  {heat(dv,x,y,train,temp,ptitle)}, 
  error=function(cond) {
  stop("Error in minimum size parameter")}
  ))

 # temproary data creation 
temp = data.frame(matrix(rep(0,length(idv)*nrow(train)),nrow(train),length(idv)))
name = names(train)
nocol = ncol(train) + 1
sub_name = idv
for (j in 1:length(sub_name)){
sub_name[j] = paste(sub_name[j],"cat",sep = "_")
}
colnames(temp) <- c(sub_name)


# Bin the categorical variables and store it in Temp
for(i in 1:length(idv)){
         if(length(levels(train[,idv[i]]))>5){
 ptitle = paste(idv[i])
 tup = xyz.bin4(train[,idv[i]],min.obs)   
temp[,i] = xyz.bin3(tup,method) 
 if (plot.bin){
 heat.map(dv,idv[i],sub_name[i],train,temp,ptitle)
 devAskNewPage(ask=TRUE)
 }
 } else{
                                                     temp = as.data.frame(temp[-i])
     warning(sprintf("Binning not done for variable %s as it has less than 5 levels",idv[i]))
 }
}
               
                train = cbind(train,temp) 
for(i in nocol:ncol(train)){
     if(length(levels(train[,i]))>32) {
 
     warning(sprintf("Even after binning variable %s has more than 32 levels",names(train)[i]))
 
 }
}
    return(train) 

 
}
