new.beta.val <- function(llold,penden.env) {
  M <- get("M",penden.env)
  N <- get("N",penden.env)
  K <- get("K",penden.env)
  Derv1.obj <- Derv1(penden.env)
  assign("Derv1.cal",Derv1.obj$Derv1.cal,penden.env)
  assign("f.hat.val",Derv1.obj$f.hat.val,penden.env)
  Derv2.obj <- Derv2(penden.env,get("lambda0",penden.env))
#####
  Derv2.invers <- -my.positive.definite.solve(Derv2.obj$Derv2.pen)
  #Derv2.pen <- Derv2.obj$Derv2.pen
  #Derv2.cal <- Derv2.obj$Derv2.cal
  beta.val <- get("beta.val",penden.env)
  direc <- Derv2.invers%*%Derv1.obj$Derv1.pen
  #val1 <- seq(1:((M-1)*N))
  #val1.2 <- seq(((M-1)*N+1),length(direc[,1]))
  #direc1 <- (direc[val1])
  #direc2 <- (direc[val1.2])
  direc.new <- c(direc[seq(1:((M-1)*N))],rep(0,N),direc[seq(((M-1)*N+1),length(direc[,1]))])
  #direc1 <- direc[seq(1:((M-1)*N))]
  #direc2 <- direc[seq(((M-1)*N+1),length(direc[,1]))]
  #direc3 <- matrix(0,1,N)
  val <- llold
  step <- 1
  beta.val <- get("beta.val",penden.env)
  #ck.temp <- ck(penden.env,(beta.val-step*c(direc1,direc3,direc2)))
  ck.temp <- ck(penden.env,(beta.val-step*direc.new))

  while(any(ck.temp=="Inf")|any(is.na(ck.temp))) {
    step <- step/2
    ck.temp <- ck(penden.env,(beta.val-step*direc.new))
  }
  f.hat.val.temp <- f.hat(penden.env,ck.temp)
  #val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val.temp,(beta.val-step*c(direc1,direc3,direc2)))
  val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val.temp,(beta.val-step*direc.new))
 
  while(val2=="Inf") {
    step <- step/2
    #ck.temp <- ck(penden.env,(beta.val-step*c(direc1,direc3,direc2)))
    ck.temp <- ck(penden.env,(beta.val-step*direc.new))
    f.hat.val.temp <- f.hat(penden.env,ck.temp)
    #val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val.temp,(beta.val-step*c(direc1,direc3,direc2)))
    val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val.temp,(beta.val-step*direc.new))
  }

  while(val2=="NaN") {
    if (step > 1e-15) {
      step <- step/2
      #ck.temp <- ck(penden.env,(beta.val-step*c(direc1,direc3,direc2)))
      ck.temp <- ck(penden.env,(beta.val-step*direc.new))
      f.hat.val.temp <- f.hat(penden.env,ck.temp)
      #val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val.temp,(beta.val-step*c(direc1,direc3,direc2)))
      val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val.temp,(beta.val-step*direc.new))
   }
    else {break}
  }
  index <- matrix(1:length(ck.temp),length(ck.temp[,1]),K,byrow=TRUE)
  obj5 <- f.hat.val.temp^2

  while(any(obj5=="Inf")) {
     step <- step/2
     ck.temp <- ck(penden.env,(beta.val-step*direc.new))
     f.hat.val.temp <- f.hat(penden.env,ck.temp)
     val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val.temp,(beta.val-step*direc.new))
     obj5 <- f.hat.val.temp
  }
  if(val2>val) calc <- FALSE else calc <- TRUE
  while(calc & step>=1e-6) {
    if((val2 <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat(penden.env,ck(penden.env,(beta.val-step*direc.new))),(beta.val-step*direc.new)))<=val) step <- step/2
    else {
      calc <- FALSE
    }
  }
  val2<- pen.log.like(penden.env,get("lambda0",penden.env),f.hat(penden.env,ck(penden.env,(beta.val-step*direc.new))))
  if(val2>val) {
    
    #beta.temp <- beta.val-step*c(direc1,direc3,direc2)
    beta.temp <- beta.val-step*direc.new
    assign("beta.val",beta.temp,penden.env)
    ck.temp <- ck(penden.env,beta.temp)
    assign("ck.temp",ck.temp,penden.env)
    f.hat.val <- f.hat(penden.env,ck.temp)
    assign("f.hat.val",f.hat(penden.env,get("ck.temp",penden.env)),penden.env)
    #assign("Derv2.cal",Derv2.cal,penden.env)
    assign("Derv2.cal",Derv2.obj$Derv2.cal,penden.env)
    #assign("Derv2.pen",Derv2.pen,penden.env)
    assign("Derv2.pen",Derv2.obj$Derv2.pen,penden.env)
    return(list(Likelie=val2,step=step))
  }
  if(val2<=val) return(list(Likelie=NA))
}
