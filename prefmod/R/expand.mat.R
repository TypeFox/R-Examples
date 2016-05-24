# rows of matrix mat are expandend according to
# frequencies given in vector freq,
#    freq must have a length of nrow(mat)

expand.mat<-function(mat,freq)
{
    mat<-as.matrix(mat)
    mat[rep(seq_len(nrow(mat)), freq), ]
}
