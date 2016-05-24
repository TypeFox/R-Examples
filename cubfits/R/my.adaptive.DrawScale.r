### Get the specific function according to the options.
get.my.update.DrawScale <- function(adaptive){
  if(!any(adaptive[1] %in% .CF.CT$adaptive)){
    stop("adaptive is not found.")
  }
  ret <- eval(parse(text = paste("my.update.DrawScale.",
                                 adaptive[1], sep = "")))
  assign("my.update.DrawScale", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.update.DrawScale().


### No adaptive.
my.update.DrawScale.none <- function(var.names, default.DrawScales){
  invisible()
} # End of my.update.DrawScale.none().


### Update scaling factors for every gene.
my.update.DrawScale.simple <- function(var.names, default.DrawScales){
  for(i in 1:length(var.names)){
    ### Update new scaling factors.
    ret <- my.DrawScale.scaling(var.names[i], .cubfitsEnv$curr.renew,
                                default.DrawScales[i])

    ### Update global.
    .cubfitsEnv$all.DrawScale[[var.names[i]]] <- ret

    ### Update curr.renew to global.
    .cubfitsEnv$DrawScale[[var.names[i]]][[.cubfitsEnv$curr.renew + 1]] <- ret
  }

  ### Update current window.
  .cubfitsEnv$curr.renew <- .cubfitsEnv$curr.renew + 1

  invisible()
} # End of my.update.DrawScale().

my.DrawScale.scaling <- function(var.name, curr.window, default.DrawScale){
  if(curr.window > 2){
    prev.scale <- .cubfitsEnv$DrawScale[[var.name]][[curr.window - 1]]
    prev.accept <- .cubfitsEnv$adaptive[[var.name]][[curr.window - 1]] /
                   .CF.AC$renew.iter
  }
  curr.scale <- .cubfitsEnv$DrawScale[[var.name]][[curr.window]]
  curr.accept <- .cubfitsEnv$adaptive[[var.name]][[curr.window]] /
                 .CF.AC$renew.iter
  ret <- curr.scale

  if(var.name == "b"){
    .cubfitsEnv$my.cat("- b: rbind(curr.scale, curr.accept)\n")
    .cubfitsEnv$my.print(rbind(curr.scale, curr.accept))
  }

  ### Smaller than the target.lower.
  id <- which(curr.accept <= .CF.AC$target.accept.lower)
  if(length(id) > 0){
    ret[id] <- curr.scale[id] * rep(.CF.AC$scale.decrease, length(id))
  }

  ### Larger than the target.upper.
  id <- which(curr.accept >= .CF.AC$target.accept.upper)
  if(length(id) > 0){
    ret[id] <- curr.scale[id] * rep(.CF.AC$scale.increase, length(id))
  }

  ### Replace too small and too large numbers, relatively.
  upper.DrawScale <- .CF.AC$sigma.upper * default.DrawScale
  lower.DrawScale <- .CF.AC$sigma.lower * default.DrawScale
  ret[ret >= upper.DrawScale] <- upper.DrawScale
  ret[ret <= lower.DrawScale] <- lower.DrawScale

  ### Check weird situations and set back to default values if applicable.
  if(curr.window > 2){
    ### Check bounds.
    id.scale.lower <- curr.scale <= lower.DrawScale &
                        prev.scale <= lower.DrawScale
    id.scale.upper <- curr.scale >= upper.DrawScale &
                        prev.scale >= upper.DrawScale
    id.accept.0 <- curr.accept == 0 & prev.accept == 0
    id.accept.1 <- curr.accept == 1 & prev.accept == 1
    accept.lower <- sum(curr.accept < .CF.AC$target.accept.lower)
    accept.upper <- sum(curr.accept > .CF.AC$target.accept.upper)

    ### Print.
    .cubfitsEnv$my.cat("- var.name: ", var.name, "\n", sep = "")
    .cubfitsEnv$my.cat("    scale bound reached #: lower = ",
                       sum(id.scale.lower), ", upper = ",
                       sum(id.scale.upper), "\n", sep = "")
    .cubfitsEnv$my.cat("    ill acceptance #: none = ",
                       sum(id.accept.0), ", all = ",
                       sum(id.accept.1), "\n", sep = "")
    .cubfitsEnv$my.cat("    acceptance NOT in range #: lower = ", accept.lower,
                       ", upper = ", accept.upper,
                       ", total = ", accept.lower + accept.upper,
                       "\n", sep = "")
  }

  ret
} # End of my.DrawScale.scaling().

