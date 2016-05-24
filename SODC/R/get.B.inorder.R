get.B.inorder <-
function(c, B, p){
        
        B.inorder=NULL
        s=NULL
        for(j in 1: p){
         s = c(s, seq(j,,p,c-1))
        }
        
        B.inorder=B[,s]
        
     return(B.inorder)
}
