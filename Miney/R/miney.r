## ################################################################################
## Program : Miney
## Language: R
## This program implements a simple version of the idea behind games
## known as "Minesweeper" on Microsoft Windows or "KMines" on KDE for
## Unix-like Operating Systems
##
## Copyright (C) 2010  Roland Rau
## 
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## ################################################################################

miney <- function(n) {
  # main function to play the game. Only parameter is 'n', a numeric
  # scalar, which indicates the size of the matrix to play, i.e. n=5
  # results in n x n = 5 x 5 matrix.
  
  create.matrix <- function(n, proportion=0.15) {
    # creates a matrix on which the game is played.
    # parameters are
    # n             numeric scalar   size of square matrix
    # proportion    numeric scalar   proportion of bombs in relation
    #                                to all fields 
    bombs <- sample(x=1:(n*n), size=floor(proportion*n*n), replace = FALSE)
    elements <- rep(0,n*n)
    elements[bombs] <- 1
    A <- matrix(elements, ncol=n, nrow=n, byrow=TRUE)
    return(A)
  }
  
  count.neighbors <- function(amatrix) {
    # given a matrix of zeros ('0') and ones ('1'), this function
    # counts the number of neighbors
    amatrix2 <- rbind(0,cbind(0,amatrix,0),0)
    resultsmatrix <- matrix(0, ncol=ncol(amatrix2), nrow=nrow(amatrix2))
    for (i in 2:(nrow(amatrix2)-1)) {
      for (j in 2:(ncol(amatrix2)-1)) {
        resultsmatrix[i,j] <- sum(amatrix2[(i-1):(i+1), (j-1):(j+1)]) - amatrix2[i,j]
      }
    }
    resultsmatrix <- resultsmatrix[-(c(1,nrow(resultsmatrix))), -(c(1,ncol(resultsmatrix)))]
    return(resultsmatrix)
  }
  
  plot.initial.matrix <- function(n, the.offset=0.05) {
    # function to plot the initial matrix, parameters are:
    # n              numeric scalar    size of matrix
    # the.offset     numeric scalar    determines how much white space
    #                                  there is between matrix
    #                                  elements.
    plot(x=1,y=1,type="n", xlab="",ylab="", axes=FALSE, xlim=c(0,n),ylim=c(0,n))
    for (i in 1:n) {
      for (j in 1:n) {
        rect(xleft=i-1+the.offset,xright=i-the.offset,
             ybottom=j-1+the.offset,ytop=j-the.offset,
             col="black")
      }
    }
  }
  
  plot.known.matrix <- function(n, the.offset=0.05) {
    # function to plot the currently known matrix, parameters are:
    # n              numeric scalar    size of matrix
    # the.offset     numeric scalar    determines how much white space
    #                                  there is between matrix
    #                                  elements.
  
    plot(x=1,y=1,type="n", xlab="",ylab="", axes=FALSE, xlim=c(0,n),ylim=c(0,n))
    for (i in 1:n) {
      for (j in 1:n) {
        if (known.matrix[i,j]==1) {
          rect(xleft=i-1+the.offset,xright=i-the.offset,
               ybottom=j-1+the.offset,ytop=j-the.offset, col="white")
          text(x=i+0.5-1, y=j+0.5-1, labels=neighbors.matrix[i,j])
        } else {
          rect(xleft=i-1+the.offset,xright=i-the.offset,
               ybottom=j-1+the.offset,ytop=j-the.offset, col="black")
        }
      }
    }
  }
  
  the.offset <- 0.05
  tick <- as.numeric(Sys.time())
  plot.initial.matrix(n=n)
  base.matrix <- create.matrix(n)
  neighbors.matrix <- count.neighbors(base.matrix)
  known.matrix <- matrix(0,ncol=n,nrow=n)
  plot.known.matrix(n=n, the.offset=the.offset)
  game.over <- FALSE
  
  while(!game.over) {
    current.click <- locator(1)
    current.x <- floor(current.click$x)+1
    current.y <- floor(current.click$y)+1
    intermediate <- known.matrix
    intermediate[current.x,current.y] <- 1
    known.matrix <- intermediate
    
    if (base.matrix[current.x,current.y]==1) {
      tock <- as.numeric(Sys.time())
      time.elapsed <- tock-tick
      rect(xleft=current.x-1+the.offset,xright=current.x-the.offset,
           ybottom=current.y-1+the.offset,ytop=current.y-the.offset, col="red")
      title(main="Game Over!\nFAILURE!!!",
            sub=paste("It took you", round(time.elapsed), "seconds!"))
      game.over <- TRUE
    } else {
      if (all(base.matrix+known.matrix==1)) {
        tock <- as.numeric(Sys.time())
        time.elapsed <- tock-tick
        plot.known.matrix(n=n, the.offset=the.offset)
        title(main="You WON!!!!!!",
              sub=paste("It took you", round(time.elapsed), "seconds!"))
        game.over <- TRUE
      } else {
        plot.known.matrix(n=n, the.offset=the.offset)
      }
    }
  }
}
