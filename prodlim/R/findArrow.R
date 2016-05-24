findArrow <- function(Box1,
                      Box2,
                      Box1Dim,
                      Box2Dim,
                      verbose=FALSE){
  left1 <- Box1[1]
  bottom1 <- Box1[2]
  left2 <- Box2[1]
  bottom2 <- Box2[2]
  width1 <- Box1Dim[1]
  height1 <- Box1Dim[2]
  width2 <- Box2Dim[1]
  height2 <- Box2Dim[2]
  
  #    ############################
  #    #p3         p4           p5#
  #    #                          #
  #    #                          #
  #    #p2                      p6# 
  #    #                          #
  #    #                          #
  #    #p1          p8          p7#
  #    ############################

  box1 <- list(left=as.numeric(left1),
               right=as.numeric(left1+width1),
               mid.horizontal=as.numeric(left1+width1/2),
               bottom=as.numeric(bottom1),
               top=as.numeric(bottom1+height1),
               mid.vertical=as.numeric(bottom1+height1/2))
  
  box1$p1 <- c(x=box1$left,y=box1$bottom)
  box1$p2 <- c(x=box1$left,y=box1$mid.vertical)
  box1$p3 <- c(x=box1$left,y=box1$top)
  box1$p4 <- c(x=box1$mid.horizontal,y=box1$top)
  box1$p5 <- c(x=box1$right,y=box1$top)
  box1$p6 <- c(x=box1$right,y=box1$mid.vertical) 
  box1$p7 <- c(x=box1$right,y=box1$bottom) 
  box1$p8 <- c(x=box1$mid.horizontal,y=box1$bottom) 
  
  box2 <- list(left=as.numeric(left2),
               right=as.numeric(left2+width2),
               mid.horizontal=as.numeric(left2+width2/2),
               bottom=as.numeric(bottom2),
               top=as.numeric(bottom2+height2),
               mid.vertical=as.numeric(bottom2+height2/2))
  box2$p1 <- c(x=box2$left,y=box2$bottom)
  box2$p2 <- c(x=box2$left,y=box2$mid.vertical)
  box2$p3 <- c(x=box2$left,y=box2$top)
  box2$p4 <- c(x=box2$mid.horizontal,y=box2$top)
  box2$p5 <- c(x=box2$right,y=box2$top)
  box2$p6 <- c(x=box2$right,y=box2$mid.vertical) 
  box2$p7 <- c(x=box2$right,y=box2$bottom) 
  box2$p8 <- c(x=box2$mid.horizontal,y=box2$bottom) 
  ##   boxwidth <- abs(box1$left-box1$right)
  ##   boxheight <- abs(box1$top-box1$bottom)
  direction <- 1
  if (box2$mid.horizontal <box1$mid.horizontal){  
    if (verbose) print("change boxes")
    direction <- 2
    tmpBox <- box1
    box1 <- box2
    box2 <- tmpBox
  }
  
  ##   points(box1$left,box1$bottom)
  ##   points(box2$left,box2$bottom,col=2)
  ##   print(list(box1,box2))
  
  if (box1$mid.horizontal==box2$mid.horizontal){
    if (box1$bottom<box2$bottom){
      #########################
      ####      2        
      ####      |
      ####      1
      ######################### 
      if (verbose==TRUE) print("case 0a: top -> bottom")
      out <- list(from=box1$p4,to=box2$p8)
    }
    else{
      #########################
      ####      1      
      ####      |
      ####      2
      ######################### 
      if (verbose==TRUE) print("case 0: bottom -> top")
      out <- list(from=box1$p8,to=box2$p4)
    }
  }
  else{
    ##     if (box1$right<=box2$left){        
    if (box1$bottom<=box2$bottom){
      if (box1$top >= box2$bottom){
        #########################
        ####          2
        ####  1   ->  
        ####
        #########################
        if (verbose==TRUE) {
          print("case 2: mid.left -> mid.right")
          print(c(from=box1$p6,to=box2$p2))
        }
        out <- list(from=box1$p6,to=box2$p2)
        #########################
        ## THIS IS A SPECIAL CASE
        ####
        ####  1   ->  2
        ####
        #########################
      }
      else{ # box1$top < box2$bottom
        if ((box2$bottom-box1$top) <= (box2$left-box1$right)){
          if ((box2$bottom-box1$top) <= .5*(box2$left-box1$right)){
            #########################
            ####    -> 2
            ####   /   
            ####  1    
            #########################
            if (verbose==TRUE) print("case 3a: corner.left.top -> mid.right")
            out <- list(from=box1$p5,to=box2$p2)
          }
          else{
            #########################
            ####     2
            ####   /   
            ####  1    
            #########################
            if (verbose==TRUE) print("case 3b: corner.left.top -> corner.right.bottom")
            out <- list(from=box1$p5,to=box2$p1)
          }
        }
        else{
          #########################
          ####    2
          ####   / 
          ####  | 
          ####  1    
          #########################
          if (verbose==TRUE) print("case 4: top.left -> bottom.right")
          out <- list(from=box1$p4,to=box2$p8)
        }
      }
    }
    ##   }
    else{ ## box1$bottom>box2$bottom
      if (box2$top>=box1$bottom){
        #########################
        ####          
        ####  1   ->  
        ####          2
        #########################
        if (verbose==TRUE) {
          print("case 5: mid.left -> mid.right")
          print(c(from=box1$p6,to=box2$p2))
        }
        
        out <- list(from=box1$p6,to=box2$p2)
      }
      else{
        if ((box1$bottom-box2$top) <= (box2$left-box1$right)){
          ## print((box1$bottom-box2$top) <= .5*(box2$left-box1$right))
          if ((box1$bottom-box2$top) <= .5*(box2$left-box1$right)){
            #########################
            ####  1        
            ####    \
            ####      2
            #########################
            if (verbose==TRUE) print("case 6a: corner.left.bottom -> mid.right")
            out <- list(from=box1$p7,to=box2$p2)
          }
          else{
            #########################
            ####  1        
            ####    \
            ####      2
            #########################
            if (verbose==TRUE) print("case 6b: corner.left.bottom -> corner.right.top")
            out <- list(from=box1$p7,to=box2$p3)
          }
        }
        else{
          if (box1$bottom>=box2$top){
            #########################
            ####  1        
            ####     \->   2
            ####
            #########################
            if (verbose==TRUE) print("case 7: top.left -> bottom.right")
            out <- list(from=box1$p8,to=box2$p4)
          }
        }
      }
    }
  }
  if (direction==2){
    names(out) <- c("to","from")
  }
  out
}
