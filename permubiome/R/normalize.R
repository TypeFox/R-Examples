normalize <-
function (numz = 0.5, method = 1)
{
load("permubiome.RData")
df_norm<-df
if (method == 1){
y <- array(, nrow(df_norm))
for (j in 1:nrow(df_norm)){
y[j]<-sum(df_norm[j,3:ncol(df_norm)])
}
for (l in 3:ncol(df_norm)){
for (m in 1:nrow(df_norm)){
df_norm[m,l]<-round((df_norm[m,l]/y[m])*1000000, digits=0)
}
}
for (i in ncol(df_norm):3) {
if (sum(df_norm[,i] == "0") >= (nrow(df_norm)*numz)) {
df_norm[,i]<-NULL
}
}
} else if (method == 2){
for (i in ncol(df_norm):3) {
if (sum(df_norm[,i] == 0) >= (nrow(df_norm)*numz)) {
df_norm[,i]<-NULL
}
}
sfactor_matrix<-matrix(, ncol = ncol(df_norm)-2, nrow = nrow(df_norm))
y <- array(, nrow(df_norm))
for (m in 1:nrow(df_norm)){
for (l in 3:ncol(df_norm)){
sfactor_matrix[m,l-2] <- signif((df_norm[m,l]/mean(df_norm[,l])), digits=3)
}
y[m]<-median(sfactor_matrix[m,1:ncol(sfactor_matrix)])
}
for (a in 3:ncol(df_norm)){
for (b in 1:nrow(df_norm)){
df_norm[b,a] <- round((df_norm[b,a]*y[b]), digits=0)
}
}
} else if (method == 3){
for (i in ncol(df_norm):3) {
if (sum(df_norm[,i] == 0) >= (nrow(df_norm)*numz)) {
df_norm[,i]<-NULL
}
}
quantil <- as.numeric(readline("Type the 'l' parameter (percentile between 0.01  and 0.99) to perform paulson's normalization (0.95 as default): "))
if (is.numeric(quantil) != TRUE & quantil > 1){
quantile <- 0.95
}
y <- array(, nrow(df_norm))
sfactor <- array(, nrow(df_norm))
for (m in 1:nrow(df_norm)){
x <- array(, ncol(df_norm)-2)
for (l in 3:ncol(df_norm)){
if (df_norm[m,l] <= quantile(df_norm[m,3:ncol(df_norm)], quantil, na.rm=T)){
x[l-2] <- df_norm[m,l]
} else {
x[l-2] <- NA
}
sfactor[m] <- sum(x, na.rm=T)
}
}
for (a in 3:ncol(df_norm)){
for (b in 1:nrow(df_norm)){
df_norm[b,a] <- round(((df_norm[b,a]/median(sfactor))*1000000), digits=0)
}
}
} else if (method == 0){
head(df_norm)
print(paste("Your dataset was not normalized according to method option: 0"))
} else { print(paste("Select and appropiate method for normalization: 1 ('proportions'), 2 ('anders'), 3('paulson'), or 0 ('none')"))
}
print(paste("Your normalized data now contains:", ncol(df_norm)-2, "normalize categories ready to analize" ))
save(df_norm, file="permubiome.RData")
}
