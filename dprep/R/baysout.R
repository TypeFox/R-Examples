baysout <-
function (D, blocks = 10, nclass=0, k = 3, num.out = 10) 
{#blocks must be greater or equal than num.out
#    require(FNN)
    if (sum(is.na(D))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    big = 1e+06
 #   print(dim(D))
    p=dim(D)[2]
 if(nclass==0)
 D = as.data.frame(D[,-p])
        else        D=D[as.numeric(factor(D[,p]))==nclass,-p ]
# print(dim(D))
    nrows = dim(D)[1]
    brows = dim(D)[1]
    c = 0
    Out = NULL
#number of blocks
    repet = ceiling(nrows/blocks)
    for (cycle in 1:repet) {
#       print(cycle)
        block.size = blocks
#       print(block.size)
        if (block.size * cycle <= brows) 
            block = (block.size * (cycle - 1) + 1):(block.size * 
                cycle)
        else {
            block = (block.size * (cycle - 1) + 1):brows
            block.size = length(block)
        }
#        print(block)
    #The block
        B = D[block, ]
    #initializing the matrix of distane of the neighbors 
        neighbors = matrix(rep(big, (block.size * k)), block.size, 
            k)
#    print(neighbors)
        rownames(neighbors) = rownames(B)
        neighbors = as.data.frame(neighbors)
            flag = 0
            reduce = 0
            removeB = rep(0, 0)
            removeN = rep(0, 0)
            for (j in 1:block.size) {
 neighbors[j,]=FNN::knnx.dist(D[-as.integer(rownames(B))[j],],B[j,],k)
# cat("c")
#                print(c)
# cat("score")
#                print(score(neighbors[j,]))
                    if (score(neighbors[j, ]) < c) {
                      removeB = cbind(removeB, j)
                      removeN = cbind(removeN, j)
                      reduce = reduce + 1
                    
                  }
            }
            if (reduce == dim(B)[1]) 
                flag = 1
            else if (reduce != 0) {
                block.size = block.size - reduce
                B = B[-removeB, ]
                neighbors = neighbors[-removeN, ]
            }
#print(neighbors)
        if (flag == 0) {
            #print(neighbors)
            Out = top(Out, neighbors, num.out)
#            print(Out)
            c = min(Out)
        }
#print(c)
    }#end of loop for blocks 
    xcoord = as.integer(rownames(Out))
    plot(Out, main = "Instances with Greatest score from K nearest neighbors", 
        ylab = "Sum of  Distances")
    text(1:num.out, Out, rownames(Out), cex = 0.6, pos = 4)
    return(Out)
}
