.lec.init <- function () {
  if(exists(".lec.Random.seed.table", envir=.GlobalEnv))
    rm(".lec.Random.seed.table", envir=.GlobalEnv)
  pos <- 1
  assign(".lec.Random.seed.table",list(Cg=matrix(0,nrow=0,ncol=6),
                                       Bg=matrix(0,nrow=0,ncol=6),
                                       Ig=matrix(0,nrow=0,ncol=6),
                                       AIP=matrix(0,nrow=0,ncol=2),
                                       name=c()), envir=as.environment(pos))
  .Call("r_create_current_stream",PACKAGE="rlecuyer")
  return(1)
}

.lec.exit <- function () {
  if(exists(".lec.Random.seed.table", envir=.GlobalEnv))
    rm(".lec.Random.seed.table", envir=.GlobalEnv)
  .Call("r_remove_current_stream",PACKAGE="rlecuyer")
  return(1)
}

.lec.which.stream <- function (name) {
  if(!exists(".lec.Random.seed.table", envir=.GlobalEnv))
    stop("No stream created yet");
  which <- which(.lec.Random.seed.table$name == name)
  if (length(which) <= 0) return (0)
  return(which[1])
}

.lec.GetStateList <- function (name) {
  which <- .lec.which.stream(name)
  if (which == 0) 
    stop(paste("No stream ",name," exists."));
  return (list(Cg=.lec.Random.seed.table$Cg[which,],Bg=.lec.Random.seed.table$Bg[which,],Ig=.lec.Random.seed.table$Ig[which,],Anti=.lec.Random.seed.table$AIP[which,1],IncPrec=.lec.Random.seed.table$AIP[which,2],name=.lec.Random.seed.table$name[which]))
}

.lec.SaveStateList <- function (stream) {
  which <- .lec.which.stream(stream$name)
  if (which == 0) 
    stop(paste("No stream ",stream$name," exists."));
  
  .lec.Random.seed.table$Cg[which,] <<- stream$Cg
  .lec.Random.seed.table$Bg[which,] <<- stream$Bg
  .lec.Random.seed.table$Ig[which,] <<- stream$Ig
  .lec.Random.seed.table$AIP[which,] <<- c(stream$Anti,stream$IncPrec)
  .lec.Random.seed.table$name[which] <<- stream$name
  return(which)
}

.lec.DoCreateStream <- function (name) {

  state <- .Call("r_create_stream", as.character(name),PACKAGE="rlecuyer")
  stream <- list(Cg=state[1:6],Bg=state[7:12],Ig=state[13:18],
                 Anti=as.integer(state[19]), IncPrec=as.integer(state[20]),
                 name=name)
  return(stream)
}

.lec.CreateStream <- function (names) {
  if(!exists(".lec.Random.seed.table", envir=.GlobalEnv))
    .lec.init()

  ln <- length(names)
  if (ln == 0)
    stop("Stream name is missing.")
  
  for (i in 1:ln) {
    name <- names[i]
  
    g <- .lec.DoCreateStream(name)
    which <- .lec.which.stream(name)
  
    if (which == 0) { # not in the table
      .lec.Random.seed.table$Cg <<- rbind(.lec.Random.seed.table$Cg,g$Cg[1:6])
      .lec.Random.seed.table$Bg <<- rbind(.lec.Random.seed.table$Bg,g$Bg)
      .lec.Random.seed.table$Ig <<- rbind(.lec.Random.seed.table$Ig,g$Ig)
      .lec.Random.seed.table$AIP <<- rbind(.lec.Random.seed.table$AIP,
                                           c(g$Anti, g$IncPrec))
      .lec.Random.seed.table$name <<- c(.lec.Random.seed.table$name, g$name)
    } else { # entry found in table
      warning("Name ", name, " already in use. Its state is overwritten.", immediate. = TRUE)
      .lec.SaveStateList(g)
    }
  }
}

.lec.WriteState <- function (names) {
  ln <- length(names)
  if (ln == 0)
    stop("Stream name is missing.")
  
  for (i in 1:ln) {
	name <- names[i]
  	g <- .lec.GetStateList(name)
  	cat("\nThe current state of the Rngstream ",name,":\n")
  	cat("Cg = ", g$Cg)
	}
}

.lec.WriteStateFull <- function (names) {
  ln <- length(names)
  if (ln == 0)
    stop("Stream name is missing.")
  
  for (i in 1:ln) {
	name <- names[i]
  	g <- .lec.GetStateList(name)
  	cat("\nThe Rngstream ",name,":\n")
  	cat("Anti = ",g$Anti)
  	cat("IncPrec = ",g$IncPrec, "\n")
  	cat("Ig = ",g$Ig, "\n")
  	cat("Bg = ",g$Bg, "\n")
  	cat("Cg = ",g$Cg, "\n")
	}
}

.lec.CheckSeed <- function(seed) {
  ll<-length(seed)
  if (ll < 6) {
    tmpseed<-12345
    seed<-c(seed,rep(tmpseed,6-ll))
  }
  if (ll > 6) {
  	seed<-seed[1:6]
    warning("Seed may not exceed the length of 6. Truncated to ", 
    		paste(seed, collapse=', '))
  }
  return(seed)
}

.lec.SetPackageSeed <- function(seed=rep(12345,6)) {
  if(!exists(".lec.Random.seed.table", envir=.GlobalEnv))
    .lec.init()
  seed <- .lec.CheckSeed(seed)
  .Call ("r_set_package_seed", as.double(seed), PACKAGE="rlecuyer")
  return(seed)
}

.lec.SetSeed <- function(name, seed=rep(12345,6)) {
  which <- .lec.which.stream(name)
  if (which == 0) 
    stop(paste("No stream ",name," exists."));
  
  seed <- .lec.CheckSeed(seed)
  g <- .lec.GetStateList(name)

  state <- .Call ("r_set_stream_seed", as.double(seed), as.double(g$Cg),
                 as.double(g$Bg),as.double(g$Ig),as.integer(g$Anti),
                 as.integer(g$IncPrec),as.character(g$name),
                 PACKAGE="rlecuyer")
  g<-list(Cg=state[1:6],Bg=state[7:12],Ig=state[13:18],
          Anti=as.integer(state[19]), IncPrec=as.integer(state[20]),
          name=name)

  .lec.SaveStateList(g)
  return(seed)
}

.lec.GetState <- function(name){
  g <- .lec.GetStateList(name)
  return(g$Cg)
}

.lec.Getuniform <- function(stream) {
  number<-.Call("r_randU01",as.double(stream$Cg),as.double(stream$Bg),
                as.double(stream$Ig),as.integer(stream$Anti),
                as.integer(stream$IncPrec),as.character(stream$name),
                PACKAGE="rlecuyer")
  return(number)
}

.lec.uniform <- function(name,n=1) {  
  g <- .lec.GetStateList(name)
  sample<-rep(0,n)
  for (i in 1:n) {
    res<-.lec.Getuniform(g)
    sample[i] <- res[[length(res)]]
    g <- list(Cg=res[1:6],Bg=res[7:12],Ig=res[13:18],
              Anti=res[19],IncPrec=res[20],name=name)
  }
  .lec.SaveStateList(g)
  return(sample)
}


.lec.uniform.int <- function (name, n=1, a=0,b=10) {
  rn <- .lec.uniform (name,n)
  rn <- as.integer(a + (b-a+1)*rn)
  return(rn)
}

.lec.DeleteStream <- function (names) {
  ln <- length(names)
  for (i in 1:ln) {
    name <- names[i]
    which <- .lec.which.stream(name)
    if (which == 0) {# not in the table
      warning("Stream ",name," does not exist.", immediate. = TRUE)
    } else { # found
      if (length(.lec.Random.seed.table$name) > 1) {
        .lec.Random.seed.table$Cg<<-.lec.Random.seed.table$Cg[-which,]
        .lec.Random.seed.table$Bg<<-.lec.Random.seed.table$Bg[-which,]
        .lec.Random.seed.table$Ig<<-.lec.Random.seed.table$Ig[-which,]
        .lec.Random.seed.table$AIP<<-.lec.Random.seed.table$AIP[-which,]
        .lec.Random.seed.table$name <<- .lec.Random.seed.table$name[-which]
      } else  {# table empty
      	pos <- 1
      	assign(".lec.Random.seed.table",list(Cg=matrix(0,nrow=0,ncol=6),
                                           Bg=matrix(0,nrow=0,ncol=6),
                                           Ig=matrix(0,nrow=0,ncol=6),
                                           AIP=matrix(0,nrow=0,ncol=2),
                                           name=c()), envir=as.environment(pos))
      }
    }
  }
}

.lec.StreamExists <- function(name) {
  which <-  which(.lec.Random.seed.table$name == name)
  if (length(which) > 0) return(TRUE)
  return(FALSE)
}

.lec.ResetStartStream <- function(name) {
  which <- .lec.which.stream(name)
  if (which == 0) 
    stop(paste("No stream ",name," exists."));
  .lec.Random.seed.table$Cg[which,] <<- .lec.Random.seed.table$Bg[which,] <<- .lec.Random.seed.table$Ig[which,]
}

.lec.ResetStartSubstream <- function(name) {
  which <- .lec.which.stream(name)
  if (which == 0) 
    stop(paste("No stream ",name," exists."));
  .lec.Random.seed.table$Cg[which,] <<- .lec.Random.seed.table$Bg[which,]
}

.lec.ResetNextSubstream <- function(name) {
  g <- .lec.GetStateList(name)  
  state <- .Call("r_reset_next_substream", as.double(g$Cg),
                 as.double(g$Bg),as.double(g$Ig),as.integer(g$Anti),
                 as.integer(g$IncPrec),as.character(g$name),
                 PACKAGE="rlecuyer")
  g<-list(Cg=state[1:6],Bg=state[7:12],Ig=state[13:18],
                 Anti=as.integer(state[19]), IncPrec=as.integer(state[20]),
                 name=name)
  .lec.SaveStateList(g)
}

.lec.SetAntithetic <- function (name, anti=FALSE) {
  which <- .lec.which.stream(name)
  if (which == 0) 
    stop(paste("No stream ",name," exists."));
  .lec.Random.seed.table$AIP[which,1] <<- as.integer(anti)
}

.lec.IncreasedPrecis <- function (name, incp=FALSE) {
  which <- .lec.which.stream(name)
  if (which == 0) 
    stop(paste("No stream ",name," exists."));
  .lec.Random.seed.table$AIP[which,2] <<- as.integer(incp)
}

.lec.AdvanceState <- function (name, e=1, c=0) {
  if (!is.numeric(e) | !is.numeric(c))
    stop(".lec.AdvanceState: Arguments must be numeric.")
  
  g <- .lec.GetStateList(name)  
  state <- .Call("r_advance_state", as.double(e), as.double(c),
                 as.double(g$Cg),
                 as.double(g$Bg),as.double(g$Ig),as.integer(g$Anti),
                 as.integer(g$IncPrec),as.character(g$name),
                 PACKAGE="rlecuyer")
  g<-list(Cg=state[1:6],Bg=state[7:12],Ig=state[13:18],
                 Anti=as.integer(state[19]), IncPrec=as.integer(state[20]),
                 name=name)
  .lec.SaveStateList(g)
  }

.lec.GetStreams <- function() {
  if(!exists(".lec.Random.seed.table", envir=.GlobalEnv))
    return(NULL)
  return(.lec.Random.seed.table$name)
}

.lec.CurrentStream <- function(name) {
  old.kind<-RNGkind()
  if (exists (".Random.seed") && floor(.Random.seed[1]/100) == 5) {
        # Leftover from have used user defined rng and forgot
        # to switch back after last run
        set.seed (0)
    }
  RNGkind ("user")
  g <- .lec.GetStateList(name)
  .Call("r_set_current_stream",as.double(g$Cg),
         as.double(g$Bg),as.double(g$Ig),as.integer(g$Anti),
         as.integer(g$IncPrec),as.character(g$name),
         PACKAGE="rlecuyer" )
  return(old.kind)
}

.lec.CurrentStreamEnd <- function(kind.old = c("Marsaglia-Multicarry",
                                     "Kinderman-Ramage")) {
   stream <- .Call("r_get_current_stream",PACKAGE="rlecuyer" )


   if (as.character(stream[[2]])=="") return(NULL) 
   g<-list(Cg=stream[[1]][1:6],Bg=stream[[1]][7:12],Ig=stream[[1]][13:18],
           Anti=as.integer(stream[[1]][19]),
           IncPrec=as.integer(stream[[1]][20]),
           name=as.character(stream[[2]]))
   .lec.SaveStateList(g)
    RNGkind (kind.old[1], kind.old[2])
   return(g$name)
}
