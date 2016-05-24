instantiate <-
function(equation, chemical_table, parentTable2, directionality, classCompound, id_col, inchi_col, smiles_col, formula_col) {
  cat(sprintf('Processing %s\n', equation))
  ## Remove pipe. Pipe indicates class compound like |Alcohols|, but chemical table doesn't have it.
  equation2 = gsub('\\|', '', equation)
  
  ## Split equation
  pattern = '([0-9]+ )(.+)'
  pattern2 = '(.+)\\(.+\\)'
  
  substrates = unlist(strsplit(equation2, directionality))[1]
  substrates = unlist(strsplit(substrates, ' \\+ '))
  
  ind_sub = grep(pattern, substrates)
  sub_coefficient = character(length(substrates))
  if(length(ind_sub) > 0) {
    sub_coefficient[ind_sub] = sub(pattern, '\\1', substrates[ind_sub])
  }
  
  ind_sub2 = grep(pattern2, substrates)
  sub_tail = character(length(substrates))
  if(length(ind_sub2) > 0) {
    sub_tail[ind_sub2] = sub(pattern2, '\\1', substrates[ind_sub2])
  }
  
  substrates = sub(pattern, '\\2', substrates)
  substrates = sub(pattern2, '\\1', substrates)
  
  products = unlist(strsplit(equation2, directionality))[2]
  products = unlist(strsplit(products, ' \\+ '))
  
  ind_pro = grep(pattern, products)
  pro_coefficient = character(length(products))
  if(length(ind_pro) > 0) {
    pro_coefficient[ind_pro] = sub(pattern, '\\1', products[ind_pro])
  }
  
  ind_pro2 = grep(pattern2, products)
  pro_tail = character(length(products))
  if(length(ind_pro2) > 0) {
    pro_tail[ind_pro2] = sub(pattern, '\\1', products[ind_pro2])
  }
  
  products = sub(pattern, '\\2', products)
  products = sub(pattern2, '\\1', products)
  
  ## Extract generic (class) compound
  class_substrate = substrates[substrates %in% classCompound]
  class_product = products[products %in% classCompound]
  
  if(length(class_substrate) == 1 & length(class_product) == 1) {
    id_col.x = paste(id_col, '.x', sep='')
    id_col.y = paste(id_col, '.y', sep='')
    
    substrates_class_coefficient = sub_coefficient[substrates == class_substrate]
    substrates_class_coefficient = as.numeric(substrates_class_coefficient)
    substrates_class_coefficient[is.na(substrates_class_coefficient)] = 1
    
    products_class_coefficient = pro_coefficient[products == class_product]
    products_class_coefficient = as.numeric(products_class_coefficient)
    products_class_coefficient[is.na(products_class_coefficient)] = 1
    
    child_substrate = .create_df(class_substrate,chemical_table, parentTable2, id_col=id_col, inchi_col=inchi_col, smiles_col)
    child_product = .create_df(class_product, chemical_table, parentTable2, id_col=id_col, inchi_col=inchi_col, smiles_col)
    
    child_substrate_c = lapply(child_substrate[[id_col]], .get.c.num, chemical_table, id_col, formula_col)
    child_product_c = lapply(child_product[[id_col]], .get.c.num, chemical_table, id_col, formula_col)
    
    child_substrate_h = lapply(child_substrate[[id_col]], .get.h.num, chemical_table, id_col, formula_col)
    child_product_h = lapply(child_product[[id_col]], .get.h.num, chemical_table, id_col, formula_col)
    
    child_substrate_o = lapply(child_substrate[[id_col]], .get.o.num, chemical_table, id_col, formula_col)
    child_product_o = lapply(child_product[[id_col]], .get.o.num, chemical_table, id_col, formula_col)
    
    child = expand.grid(child_substrate[[id_col]], child_product[[id_col]], stringsAsFactors=FALSE)
    child_c = expand.grid(substrates_class_coefficient * unlist(child_substrate_c), products_class_coefficient * unlist(child_product_c))
    child_h = expand.grid(substrates_class_coefficient * unlist(child_substrate_h), products_class_coefficient * unlist(child_product_h))
    child_o = expand.grid(substrates_class_coefficient * unlist(child_substrate_o), products_class_coefficient * unlist(child_product_o))
    
    colnames(child) = c(id_col.x, id_col.y)
    colnames(child_c) = c(id_col.x, id_col.y)
    colnames(child_h) = c(id_col.x, id_col.y)
    colnames(child_o) = c(id_col.x, id_col.y)
    
    if(nrow(child) == 0) {
      return(0)
    }
    
    # Fiter with InChI if they have same number of carbons
    substrates_part = substrates[substrates != class_substrate]
    products_part = products[products != class_product]
    
    substrates_part_coefficient = sub_coefficient[substrates != class_substrate]
    products_part_coefficient = pro_coefficient[products != class_product]
    
    substrates_part_coefficient = as.numeric(substrates_part_coefficient)
    substrates_part_coefficient[is.na(substrates_part_coefficient)] = 1
    
    products_part_coefficient = as.numeric(products_part_coefficient)
    products_part_coefficient[is.na(products_part_coefficient)] = 1
    
    substrates_part2 = rep(substrates_part, substrates_part_coefficient)
    products_part2 = rep(products_part, products_part_coefficient)
    
    substrates_c = lapply(substrates_part2, .get.c.num, chemical_table, id_col, formula_col)
    substrates_c = sum(unlist(substrates_c))
    
    products_c = lapply(products_part2, .get.c.num, chemical_table, id_col, formula_col)
    products_c = sum(unlist(products_c))
    
    child2 = c()
    if(substrates_c == products_c) {
      child2 = merge(child_substrate, child_product, by=inchi_col, incomparables='')
    } 
    
    if(nrow(child2) == 0) {
      # Filter with the number of C
      ind_c_balanced = which(child_c[[id_col.x]] + substrates_c == child_c[[id_col.y]] + products_c)
      
      # Filter with the number of H
      substrates_h = lapply(substrates_part2, .get.h.num, chemical_table, id_col, formula_col)
      substrates_h = sum(unlist(substrates_h))
      
      products_h = lapply(products_part2, .get.h.num, chemical_table, id_col, formula_col)
      products_h = sum(unlist(products_h))
      
      ind_h_balanced = which(child_h[[id_col.x]] + substrates_h == child_h[[id_col.y]] + products_h)
      
      # Filter with the number of O
      substrates_o = lapply(substrates_part2, .get.o.num, chemical_table, id_col, formula_col)
      substrates_o = sum(unlist(substrates_o))
      
      products_o = lapply(products_part2, .get.o.num, chemical_table, id_col, formula_col)
      products_o = sum(unlist(products_o))
      
      ind_o_balanced = which(child_o[[id_col.x]] + substrates_o == child_o[[id_col.y]] + products_o)
      
      ind_balanced = Reduce(intersect, list(ind_c_balanced, ind_h_balanced, ind_o_balanced))
      child2 = child[ind_balanced, ]
    }
    ind_same = which(child2[[id_col.x]] == child2[[id_col.y]])
    
    if(length(ind_same) > 0) {
      child3 = child2[-ind_same,]
    } else {
      child3 = child2	
    }
    
    if(nrow(child3) == 0) {
      cat(sprintf('%s cannot be instantiated because there might be lack of proper pair of children\n', equation))
      return(0)
    }
    
    # Assemble equation
    substrates_coeff = sub_coefficient[substrates != class_substrate]
    substrates_tail = sub_tail[substrates != class_substrate]
    
    substrates_part = paste(substrates_coeff, substrates_part, substrates_tail, sep='')
    
    if(length(substrates_part) > 1) {
      substrates_part = paste(substrates_part, collapse=' + ')
    }
    instance_coeff = sub_coefficient[substrates == class_substrate]
    instance_tail = sub_tail[substrates == class_substrate]
    
    child_tmp = child3
    child_tmp[[id_col.x]] = paste(instance_coeff, child_tmp[[id_col.x]], instance_tail, sep='')
    
    sub_instance_table = cbind(substrates_part, child_tmp[[id_col.x]])
    sub_equation = apply(sub_instance_table, 1, paste, collapse=' + ')
    
    products_coeff = pro_coefficient[products != class_product]
    products_tail = pro_tail[products != class_product]
    
    products_part = paste(products_coeff, products_part, products_tail, sep='')
    
    if(length(products_part) > 1) {
      products_part = paste(products_part, collapse=' + ')
    }
    instance_coeff = pro_coefficient[products == class_product]
    instance_tail = pro_tail[products == class_product]
    child_tmp[[id_col.y]] = paste(instance_coeff, child_tmp[[id_col.y]], instance_tail, sep='')
    
    pro_instance_table = cbind(products_part, child_tmp[[id_col.y]])
    pro_equation = apply(pro_instance_table, 1, paste, collapse=' + ')
    instance_equation = apply(cbind(sub_equation, pro_equation), 1, paste, collapse=directionality)
    ind_balanced = check.mass.balance(instance_equation, chemical_table, id_col, formula_col, direction_type=directionality)
    ind_balanced = grep('TRUE', ind_balanced)
    child4 = child3[ind_balanced,]
    # One instance substrate has only one instance product, and it is should be most similar chemical pair
    unique_list.x = names(which(table(child4[[id_col.x]]) == 1))
    unique_table.x = child4[child4[[id_col.x]] %in% unique_list.x,c(id_col.x, id_col.y)]
    duplication_list.x = names(which(table(child4[[id_col.x]]) > 1))
    for(k in duplication_list.x) {
      target = child4[child4[[id_col.x]] == k, id_col.y]
      similarity = numeric(length(target))
      for(j in 1:length(target)) {
        sdf1 = smiles2sdf(chemical_table[chemical_table[[id_col]] == target[j],smiles_col])
        sdf2 = smiles2sdf(chemical_table[chemical_table[[id_col]] == k, smiles_col])
        similarity[j] = fmcs(sdf1, sdf2, fast=T)[['Overlap_Coefficient']]
      }
      max_target = which(similarity == max(similarity))
      if(length(max_target) > 1) {
        cat(sprintf('WARNING: For compound %s, target %s have same similarity, so you have to check manually\n', equation, paste(target[max_target], collapse=' and ')))
      }
      instance_pair = cbind(k, target[max_target])
      colnames(instance_pair) = colnames(unique_table.x)
      unique_table.x = rbind(unique_table.x, instance_pair)
    }
    
    unique_list.y = names(which(table(unique_table.x[[id_col.y]]) == 1))
    unique_table = unique_table.x[unique_table.x[[id_col.y]] %in% unique_list.y,c(id_col.x,id_col.y)]
    if(nrow(unique_table)>0) {
      for(k in 1:nrow(unique_table)) {
        tmpChild1 = unique_table[k,id_col.x]
        tmpChild2 = unique_table[k,id_col.y]
        sdf1 = smiles2sdf(chemical_table[chemical_table[[id_col]] == tmpChild1,smiles_col])
        sdf2 = smiles2sdf(chemical_table[chemical_table[[id_col]] == tmpChild2, smiles_col])
        ap1 = sdf2ap(sdf1)
        cid(ap1) = 'tmp1'
        ap2 = sdf2ap(sdf2)
        cid(ap2) = 'tmp2'
        fp1 = desc2fp(ap1, descnames=1024, type="FPset")
        fp2 = desc2fp(ap2, descnames=1024, type="FPset")
        tmpSimilarity = fpSim(x=fp1, y=fp2, method='Tanimoto')
        if(tmpSimilarity < 0.5) {
          cat(sprintf('WARNING: %s and %s have similarity %s, you should check manually\n', tmpChild1, tmpChild2, tmpSimilarity))
        }
      }
    }
    
    duplication_list.y = names(which(table(unique_table.x[[id_col.y]]) > 1))
    for(k in duplication_list.y) {
      target = unique_table.x[unique_table.x[[id_col.y]] == k, id_col.x]
      similarity = numeric(length(target))
      for(j in 1:length(target)) {
        sdf1 = smiles2sdf(chemical_table[chemical_table[[id_col]] == target[j],smiles_col])
        sdf2 = smiles2sdf(chemical_table[chemical_table[[id_col]] == k, smiles_col])
        similarity[j] = fmcs(sdf1, sdf2, fast=T)[['Overlap_Coefficient']]
      }
      max_target = which(similarity == max(similarity))
      if(length(max_target) > 1) {
        message(sprintf('WARNING: For compound %s, target %s have same similarity, so you have to check manually', k, paste(target[max_target], collapse=' and ')))
      }
      if(max(similarity) < 0.5) {
        message(sprintf('WARNING: %s and %s have similarity %s, you should check manually', k, target[j], similarity[similarity == max(similarity)]))
      }
      instance_pair = expand.grid(target[max_target],k)
      colnames(instance_pair) = c(id_col.x, id_col.y)
      unique_table = rbind(unique_table, instance_pair)
    }
    
    ind_final = logical(nrow(child4))
    for(l in 1:nrow(child4)) {
      index_x = which(child4[l,id_col.x] == unique_table[id_col.x])
      if(length(index_x) > 0 && child4[l, id_col.y] == unique_table[index_x, id_col.y]) {
        ind_final[l] = TRUE
      }
    }
    
    instance_equation_final = instance_equation[ind_balanced][ind_final]
    return(instance_equation_final)
  } else {
    message(sprintf('Equation: %s is not instantiated by this function', equation))
    return(0)
  }
}
