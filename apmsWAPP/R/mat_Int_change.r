##### Translation Count-Matrix <-> Interactiontable (format required for SAINT)

####     Function Interactiontable -> Matrix:
# Input: Interactiontable as required in SAINT (including zero counts)

int2mat <- function(IntSaint) {
prots <- as.character(unique(IntSaint$V3))
                                              #transfer inttab -> inttab.mat (matrix)
rowspec <- function(i){                       # = matrix row i
protframe <-IntSaint[grep(paste("^",i,"$",sep=""), IntSaint$V3), ]
return(protframe$V4[match(as.character(unique(IntSaint$V1)),protframe$V1)])}
# spectral counts in matrix-form:
IntSaint.mat <- t(vapply(prots, rowspec , rep(0, times=length(as.character(unique(IntSaint$V1))))))
colnames(IntSaint.mat) <- as.character(unique(IntSaint$V1))
return(IntSaint.mat)
}

# example
# baitint.mat <- int2mat("bait_LBcl_IntSaint.txt")


#####     Function Matrix -> Interactiontable
# Input: matrix count table, baittable as classifier for ctrls and bait-samples

mat2int <- function(mat, baittab){
inttab <- data.frame(V1=NA, V2=NA, V3=NA, V4=0.0)
for ( i in 1:dim(mat)[2]) {
V1 <- rep(colnames(mat)[i], times= dim(mat)[1] )
V2 <- rep(as.character(baittab$V2[match(colnames(mat)[i], baittab$V1)]), times=dim(mat)[1] )
V3 <- rownames(mat)
V4 <- mat[,i]
inttab <- rbind(inttab, cbind(V1,V2,V3,V4))
}
inttab <- inttab[-1,]
inttab$V4 <- as.numeric(inttab$V4)
rownames(inttab) <- NULL
return(inttab)
}

