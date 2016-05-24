morph <- function(data, frames = 25, modnames = 0, Xrange = 0, Yrange = 0) {
  ### This function produced a morphing from a starting point to a finishing point
  ### 'data' is for coordinates, X1, X2, Y1 e Y2
  
  ### 'frames' indicates the number of intermediate positions of the morphing (the smoothness)
  n = nrow(data)
  
  if (Xrange[1] == 0) {
    mins = apply(data, 2, min)
    maxs = apply(data, 2, max)
    
    minX = min(mins[1:2])
    maxX = max(maxs[1:2])
    Xrange = c(minX, maxX)
    minY = min(mins[3:4])
    maxY = max(maxs[3:4])
    Yrange = c(minY, maxY)
  }
  #plot(data[1,1:2],data[1,3:4],xlim=c(minX,maxX),ylim=c(minY,maxY),type='n')
  
  out = list()
  out$Xs = matrix(0, nrow = frames, n)
  out$Ys = matrix(0, nrow = frames, n)
  
  X = matrix(c(0), nrow = n, ncol = 2)
  Y = matrix(c(0), nrow = n, ncol = 2)
  
  for (i in 1:n) {
    
    
    X[i, ] = data[i, 1:2]
    Y[i, ] = data[i, 3:4]
    
    if ((X[i, 1] != X[i, 2]) | (Y[i, 1] != Y[i, 2])) {
      #mormod=lm(c(Y[i,y1],Y[i,y2])~c(x1,x2))
      mormod = lm(Y[i, ] ~ X[i, ])
      # print(mormod$coefficients)
      bzero = mormod$coefficients[1]
      b_one = mormod$coefficients[2]
      
      xfra = seq(from = X[i, 1], to = X[i, 2], length.out = frames)
      yfra = bzero + b_one * xfra
      
    }
    if ((X[i, 1] == X[i, 2]) & (Y[i, 1] == Y[i, 2])) {
      
      xfra = matrix(X[i, 1], nrow = 1, ncol = frames)
      yfra = matrix(Y[i, 1], nrow = 1, ncol = frames)
    }
    out$Xs[, i] = xfra
    out$Ys[, i] = yfra
    
    
  }
  if (modnames[1] == 0) {
    modnames = LETTERS[1:n]
  }
  out
} 