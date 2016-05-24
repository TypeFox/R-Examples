shift.mean <-
function(predator.frame, prey.frame) {
  ## first check if the system is feasible. if it's not, return -1
  d = dim(predator.frame)[2]
  if (d == 1) {
    ## in this case the output is easy
    return(colMeans(prey.frame) - predator.frame[1,])
  }
  tmp = isotopesolve(predator.frame, prey.frame)
  retval = 0
  if (tmp$status != 2) {
    out = construct.linearsys(predator.frame, prey.frame)
    ## produce the vertices of the above linear system
    #print("Input to vertexenum and Vertices from C:")
    vertices = enumerate.vertices(out$A, out$b)
    #print("Vertices from R:")
    #print(vertices)
    num.vertices = dim(vertices)[1]
    ## make a bunch of random draws from a uniform dirichlet distribution
    ## make shift vectors of them by representing as a sum of vertices
    retval = colMeans(vertices)
  } else {
    retval = c()
  }
  retval
}
