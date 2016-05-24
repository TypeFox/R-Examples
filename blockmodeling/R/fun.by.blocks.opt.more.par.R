"fun.by.blocks.opt.more.par" <-
function(
  x,	#an object of class "opt.more.par"
  which=1,	#which best solution/partition should be used
  ...	#aditional parameters to function "fun.by.blocks"
){
	if(which>length(x$best)){
		which<-1
		warning("Only",length(x$best),"solutions exists. The first solution will be used.")
	}
	fun.by.blocks(M=x$M, clu=clu(x,which=which),...)
}

