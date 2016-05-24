InitXi <-
function(X,bounds=c(-0.5,0.29,0.79)){
s=dim(X)[1]
m=dim(X)[2]
initial.xi=matrix(1,s,m)
initial.xi[which(X>bounds[1])]=2
initial.xi[which(X>bounds[2])]=3
initial.xi[which(X>bounds[3])]=4
return(initial.xi)
}
