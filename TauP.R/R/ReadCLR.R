ReadCLR <-
function(filename, z = 'default'){
  model = EmptyModel()
  layers = list()
  S = scan(filename, what = character(0), comment.char = '/', sep = '\n')
  n = length(S)
  j = 0
  for(i in 1:n){
    line = strsplit(S[i],' +')[[1]]
    k = length(line)
    
    if(k == 0){
      next
    }
    
    if(line[1] == '!name'){
      model$name = paste(line[2:k], collapse = ' ')
    }else if(line[1] == '!year'){
      model$year = as.numeric(line[2])
      
    }else if(line[1] == '!planet'){
      if(line[2] == '!name'){
        model$planetname = paste(line[3:k], collapse = ' ')
      }else if(line[2] == '!radius'){
        model$rp = as.numeric(line[3])
      }
      
    }else if(line[1] == '!layer'){
      if(line[2] == '!start'){
        temp.layer = list()
        temp.layer$name = paste(line[3:k], collapse = ' ')
        temp.layer$vp = NaN
        temp.layer$vs = NaN
        temp.layer$rho = NaN
        temp.layer$qp = NaN
        temp.layer$qs = NaN
      }else if(line[2] == '!depth'){
        temp.layer$depth = as.numeric(line[3:4])
      }else if(line[2] == '!vp'){
        temp.layer$vp = as.numeric(line[3:k])
      }else if(line[2] == '!vs'){
        temp.layer$vs = as.numeric(line[3:k])
      }else if(line[2] == '!rho'){
        temp.layer$rho = as.numeric(line[3:k])
      }else if(line[2] == '!qp'){
        temp.layer$qp = as.numeric(line[3:k])
      }else if(line[2] == '!qs'){
        temp.layer$qs = as.numeric(line[3:k])
      }else if(line[2] == '!end'){
        j = j + 1
        layers[[j]] = temp.layer
      }
    
    }else if(line[1] =='!discon'){
      if(line[2] == '!depth'){
        depth = as.numeric(line[3])
        name = tolower(paste(line[4:k], collapse = ' '))
        if(grepl('conr', name)){
          model$conr = depth
        }else if(grepl('moho', name) | grepl('upper.mantle', name)){
          model$moho = depth
        }else if(grepl('410', name) | grepl('olivine.alpha.beta', name) | grepl('transition.zone', name)){
          model$d410 = depth
        }else if(grepl('520', name) | grepl('olivine.beta.gamma', name)){
          model$d520 = depth
        }else if(grepl('660', name) | grepl('670', name) | grepl('olivine.gamma.perovskite', name) | grepl('lower.mantle', name)){
          model$d660 = depth
        }else if(grepl('cmb', name) | grepl('core.mantle', name) | grepl('mantle.core', name) | grepl('outer.core', name)){
          model$cmb = depth
        }else if(grepl('icb', name) | grepl('inner.core', name)){
          model$icb = depth
        }
      }
    } 
  
  } 
  n = length(layers)
  if(z == 'default'){
    z = seq(from = 0, to = model$rp, by = 20)
  }
  
  for(i in 1:n){
    layer = layers[[i]]
    layer$z = sort(unique(c(layer$depth, z[z < max(layer$depth) & z > min(layer$depth)])))
    rn = 1 - layer$z/model$rp 
    vp = 0
    vs = 0
    rho = 0
    qp = 0
    qs = 0
    for(j in 1:length(layer$vp)){
      vp = vp + rn^(j-1) * layer$vp[j]
    }
    for(j in 1:length(layer$vs)){
      vs = vs + rn^(j-1) * layer$vs[j]
    }
    for(j in 1:length(layer$rho)){
      rho = rho + rn^(j-1) * layer$rho[j]
    }
    for(j in 1:length(layer$qp)){
      qp = qp + rn^(j-1) * layer$qp[j]
    }
    for(j in 1:length(layer$qs)){
      qs = qs + rn^(j-1) * layer$qs[j]
    }
    model$z = c(model$z, layer$z)
    model$vp = c(model$vp, vp)
    model$vs = c(model$vs, vs)
    model$rho = c(model$rho, rho)
    model$qp = c(model$qp, qp)
    model$qs = c(model$qs, qs)
   }
  
  return(model)
}

