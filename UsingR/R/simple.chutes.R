##' simulate chutes and ladders
##'
##' @param sim logical
##' @param return.cl return values
##' @param cl make the game
##'
##' @export
"simple.chutes" <-
  function(sim=FALSE,return.cl=FALSE,cl=make.cl()) {
    ## make a chutes and ladder board (return.cl) or
    ## simulate a trajectory (sim=TRUE)
    ## board from http://www.ahs.uwaterloo.ca/~musuem/vexhibit/Whitehill/snakes/snakes.gif
    ## to create multiple ones, pass in cl for less overhead. eg.
    ## cl<-simple.chutes(return.cl=TRUE)
    ## results<-c();for(i in 1:1000) results[i]<-simple.chutes(sim=TRUE,cl)

    make.cl <- function() {
      ## make the cl matrix. Needs the board.
      cl<-matrix(0,nrow=100,ncol=100)
      for(i in 1:94) cl[i,i+1:6]<-rep(1/6,6)
      for(i in 95:100) cl[i,i:100]<-c((i-94)*1/6,rep(1/6,100-i))
      ## now for chutes and ladders
      from.to<-matrix(c(
                        44,1,
                        37,6,
                        32,10,
                        15,71,
                        22,58,
                        89,46,
                        99,56,
                        95,65,
                        66,96),nrow=2)
      for(i in 1:dim(from.to)[2]) {
        from<-from.to[1,i];to<-from.to[2,i]
        cl[,to]<-cl[,to]+cl[,from]
        cl[,from]<-0*cl[,from]
      }
      cl
    }

    ## lazy eval. Pass in for speed
    sim.cl <- function(cl) {
      traj<-c(1)                        #start at 1
      n<-1
      while(traj[n] != 100) {
        traj[n+1]<-sample(1:100,1,prob=cl[traj[n],])
        n<-n+1
      }
      traj                          
    }

    ## now what to do?
    if (sim) return(sim.cl(cl))
    if (return.cl) make.cl()
    
  }
