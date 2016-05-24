plot_beta <- function( beta, bpmap=NA, block.bounds.vec = NA ){
  if(!is.na(bpmap)){
    plot( bpmap, ifelse(abs(beta)>1e-6,abs(beta),0), pch = 20, col= 1, type = "p", ylim = c(0,max(abs(beta))), xlab = "kB", ylab = "|beta|" )
  }else{
    plot(1:length(beta), ifelse(abs(beta)>1e-6,abs(beta),0), pch = 20, col= 1, type = "p", ylim = c(0,max(abs(beta))), xlab = "SNP", ylab = "|beta|" )
  }
  if(!is.na(block.bounds.vec) & !is.na(bpmap) ){
    for( block.bound in block.bounds.vec ){
      lines( rep( block.bound, 2), c(0, max(abs(beta)) ), lty = 3 )
    }
  }
}
