
# 
# box_product_sparse<-function(B1, B2){
#   # troubleshoot the input matrices
#   classVector <- c("spam", "Matrix", "matrix") 
#   if((class(B1) == "spam")|(class(B2) == "spam")){
#     if((class(B2) == "Matrix")|(class(B2) == "Matrix")){
#       stop("Input matrices have differing sparse formats")
#     } else {outClass <- "spam"}
#   }
#   if((class(B1) == "Matrix")|(class(B2) == "Matrix")) outClass <- "Matrix"
#   if((class(B1) == "matrix") & (class(B2) == "matrix")) outClass <- "spam"
#   if((!class(B1) %in% classVector)|(!class(B2) %in% classVector)){
#     stop("Input objects not of a recognised matrix format")
#   }
#     
# 
#   n<-nrow(B1)
#   n.a<-ncol(B1)
#   n.b<-ncol(B2)
#   
#   if(class(B1) == "Matrix"){
#     a<-which(B1>0, arr.ind = T)
#     ind<-order(a[,1])
#     a<-a[ind,]
#     a.row<-a[,1]
#     a.col<-a[,2]
#     a.val<-B1[a]
#     by1<-rowSums(B1>0)
#   }
#   if(class(B1) == "spam"){
#     a<-triplet(B1)
#     a.row<-a$indices[,1]
#     a.col<-a$indices[,2]
#     a.val<-a$values
#     by1<-as.vector(summary(as.factor(a.row)))
#   }
#   if(class(B2) == "Matrix"){
#     b<-which(B2>0, arr.ind = T)
#     ind<-order(b[,1])
#     b<-b[ind,]
#     b.row<-b[,1]
#     b.col<-b[,2]
#     b.val<-B2[b]
#     by2<-rowSums(B2>0)
#   }
#   if(class(B2) == "spam"){
#     b<-triplet(B2)
#     b.row<-b$indices[,1]
#     b.col<-b$indices[,2]
#     b.val<-b$values
#     by2<-as.vector(summary(as.factor(b.row)))
#   }
#   
#   prod.row<-rep(a.row, rep(by2, by1))
#   num.list<-split(b.col, b.row)
#   for(i in 1:length(num.list)){num.list[[i]]<-rep(num.list[[i]], by1[i])}
#   prod.col<-(rep(a.col, rep(by2, by1))-1)*(n.b) + unlist(num.list)
#   num.list<-split(b.val, b.row)
#   for(i in 1:length(num.list)){num.list[[i]]<-rep(num.list[[i]], by1[i])}
#   prod.val<-rep(a.val, rep(by2, by1))*unlist(num.list)
#   if(outClass == "spam") out <- spam(list(i = prod.row, j = prod.col, prod.val), nrow = n, ncol = n.a*n.b)   
#   if(outClass == "Matrix") out <- sparseMatrix(i=prod.row, j = prod.col, x = prod.val, dims = c(n, n.a*n.b)) 
#   out
# }

not_sparse_box_product<-function(A, B){
  unit.vec<-matrix(1, nrow = 1, ncol = ncol(B))
  (A %x% unit.vec) * (unit.vec %x% B)}
