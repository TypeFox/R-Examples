instantiate.rhea = function(equation, chemical_table, id_col = 'chebi', parent_col='parent', formula_col='formula', smiles_col='smiles', inchi_col='inchi', direction_type=c(' <\\?> ', ' <=> ', ' => ')) {
  result = list()
  no_pair = 0
  # Substitute 1 for n
  pattern_localization = '\\(.+\\)'
  testEquation = gsub(pattern_localization, '', equation)
  
  ind_n = grep('n', testEquation)
  if(length(ind_n) >0) {
    message(sprintf('%s coefficient "n" containing equations are found being processed by substituting 1-10 for n...', length(ind_n)))
    for(n in 1:10) {
      tmpEquation = gsub('\\(n\\+2\\)', n+2, testEquation[ind_n])
      tmpEquation = gsub('\\(n\\+1\\)', n+1, tmpEquation[ind_n])
      tmpEquation = gsub('7n', 7*n, tmpEquation)
      tmpEquation = gsub('2n', 2*n, tmpEquation)
      tmpEquation = gsub('n ', paste(n, ' ', sep=''), tmpEquation)
      equation = c(equation, tmpEquation)
    }
    equation = equation[-ind_n]
  }
  # Get Participant
  participants = .get.participant(equation)
  
  # Get compound class (in this version, only considered alkyl and partly polymerized compound)
  alkyl = chemical_table[chemical_table[[id_col]] %in% participants,][grep('R', chemical_table[chemical_table[[id_col]] %in% participants, formula_col]),id_col]
  partly_polymerized = chemical_table[chemical_table[[id_col]] %in% participants,][grep('\\)n', chemical_table[chemical_table[[id_col]] %in% participants, formula_col]),id_col]
  classCompound = unique(c(alkyl, partly_polymerized))
  
  if(length(classCompound) == 0) {
    message('There is no generic compound in your equations')
    result_list = list()
    for(eq in equation) {
      result_list[[eq]] = 0
    }
    return(result_list)
  }
  cat(sprintf('%i alkyl compounds, %i partly polymerized compounds found\n', length(alkyl), length(partly_polymerized)))
  
  # Build parent-child table
  cat('Build parent table\n')
  parentTable = build.subtable(chemical_table, id_col, parent_col, '///')
  parentTable[is.na(parentTable)] = ''
  colnames(parentTable) = c("child", 'parent')
  
  # Remove stereo-related parent-child
  tmpChild = data.frame(parentTable[['child']], stringsAsFactors=F)
  colnames(tmpChild) = id_col
  
  tmpParent = data.frame(parentTable[['parent']], stringsAsFactors=F)
  colnames(tmpParent) = id_col
  
  child_smiles = join(tmpChild, chemical_table[,c(id_col, smiles_col)], by=id_col)
  child_smiles[is.na(child_smiles)] = ''
  parent_smiles = join(tmpParent, chemical_table[,c(id_col, smiles_col)], by=id_col)
  parent_smiles[is.na(parent_smiles)] = ''
  
  ind_no_smiles = grep('^$', parent_smiles[[smiles_col]])
  ind_star = grep('\\*', parent_smiles[[smiles_col]])
  
  ind_no_smiles_or_start = c(ind_no_smiles, ind_star)
  
  parentTable_1 = parentTable[ind_no_smiles_or_start,]
  child_number_of_at = str_count(child_smiles[-ind_no_smiles_or_start, smiles_col], '@')
  parent_number_of_at = str_count(parent_smiles[-ind_no_smiles_or_start, smiles_col], '@')
  
  parentTable_2 = parentTable[-ind_no_smiles_or_start,][child_number_of_at <= parent_number_of_at,]
  parentTable = rbind(parentTable_1, parentTable_2)
  
  parent_participant = participants[participants %in% unique(parentTable$parent)]
  parent_list = classCompound[classCompound %in% parent_participant]
  
  # Instantiation
  ## Determine direction symbol
  directionalities = numeric(length(equation))
  for(i in direction_type) {
    ind_direction = grepl(i, equation) 
    directionalities[ind_direction] = i
  }
  
  number_of_instantiated_generic = c()
  for(i in 1:length(equation)) {
    cat(sprintf('Processing %s\n', equation[i]))
    ## Split equation
    pattern = '([0-9]+ )(.+)'
    pattern2 = '(.+)(\\(.+\\))'
    
    substrates = unlist(strsplit(equation[i], directionalities[i]))[1]
    substrates = unlist(strsplit(substrates, ' \\+ '))
    
    ind_sub = grep(pattern, substrates)
    sub_coefficient = character(length(substrates))
    if(length(ind_sub) > 0) {
      sub_coefficient[ind_sub] = sub(pattern, '\\1', substrates[ind_sub])
    }
    
    ind_sub2 = grep(pattern2, substrates)
    sub_tail = character(length(substrates))
    if(length(ind_sub2) > 0) {
      sub_tail[ind_sub2] = sub(pattern2, '\\2', substrates[ind_sub2])
    }
    
    substrates = sub(pattern, '\\2', substrates)
    substrates = sub(pattern2, '\\1', substrates)
    
    products = unlist(strsplit(equation[i], directionalities[i]))[2]
    products = unlist(strsplit(products, ' \\+ '))
    
    ind_pro = grep(pattern, products)
    pro_coefficient = character(length(products))
    if(length(ind_pro) > 0) {
      pro_coefficient[ind_pro] = sub(pattern, '\\1', products[ind_pro])
    }
    
    ind_pro2 = grep(pattern2, products)
    pro_tail = character(length(products))
    if(length(ind_pro2) > 0) {
      pro_tail[ind_pro2] = sub(pattern2, '\\2', products[ind_pro2])
    }
    
    products = sub(pattern, '\\2', products)
    products = sub(pattern2, '\\1', products)
    
    ## Extract generic (class) compound
    class_substrate = substrates[substrates %in% parent_list]
    class_product = products[products %in% parent_list]
    if(length(class_substrate) == 1 & length(class_product) == 1) {
      number_of_instantiated_generic = c(number_of_instantiated_generic, equation[i])
      id_col.x = paste(id_col, '.x', sep='')
      id_col.y = paste(id_col, '.y', sep='')
      
      substrates_class_coefficient = sub_coefficient[substrates == class_substrate]
      substrates_class_coefficient = as.numeric(substrates_class_coefficient)
      substrates_class_coefficient[is.na(substrates_class_coefficient)] = 1
      
      products_class_coefficient = pro_coefficient[products == class_product]
      products_class_coefficient = as.numeric(products_class_coefficient)
      products_class_coefficient[is.na(products_class_coefficient)] = 1
      
      child_substrate = .create_df(class_substrate,chemical_table, parentTable, id_col=id_col, inchi_col=inchi_col, smiles_col)
      child_product = .create_df(class_product, chemical_table, parentTable, id_col=id_col, inchi_col=inchi_col, smiles_col)
      
      child_substrate_c = lapply(child_substrate[[id_col]], .get.c.num, chemical_table, id_col, formula_col)
      child_product_c = lapply(child_product[[id_col]], .get.c.num, chemical_table, id_col, formula_col)
      
      child_substrate_h = lapply(child_substrate[[id_col]], .get.h.num, chemical_table, id_col, formula_col)
      child_product_h = lapply(child_product[[id_col]], .get.h.num, chemical_table, id_col, formula_col)
      
      child_substrate_o = lapply(child_substrate[[id_col]], .get.o.num, chemical_table, id_col, formula_col)
      child_product_o = lapply(child_product[[id_col]], .get.o.num, chemical_table, id_col, formula_col)
      
      child = expand.grid(child_substrate[[id_col]], child_product[[id_col]])
      child_c = expand.grid(substrates_class_coefficient * unlist(child_substrate_c), products_class_coefficient * unlist(child_product_c))
      child_h = expand.grid(substrates_class_coefficient * unlist(child_substrate_h), products_class_coefficient * unlist(child_product_h))
      child_o = expand.grid(substrates_class_coefficient * unlist(child_substrate_o), products_class_coefficient * unlist(child_product_o))
      
      colnames(child) = c(id_col.x, id_col.y)
      colnames(child_c) = c(id_col.x, id_col.y)
      colnames(child_h) = c(id_col.x, id_col.y)
      colnames(child_o) = c(id_col.x, id_col.y)
      
      if(nrow(child) == 0) {
        result[[equation[i]]] = 0
        next    
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
      
      if(substrates_c == products_c) {
        child2 = merge(child_substrate, child_product, by=inchi_col, incomparables='')
      } else {
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
      
      if(length(ind_same) > 0 && class_substrate != class_product) {
        child3 = child2[-ind_same,]
      } else {
        child3 = child2 
      }
      if(nrow(child3) == 0) {
        message(sprintf('WARNING: %s cannot be instantiated because there might be lack of proper pair of children', equation[i]))
        result[[equation[i]]] = 0
        no_pair = no_pair + 1
        next
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
      directionalities[i] = gsub('\\\\','',directionalities[i])
      instance_equation = apply(cbind(sub_equation, pro_equation), 1, paste, collapse=directionalities[i])
      ind_balanced = check.mass.balance(instance_equation, chemical_table, id_col, formula_col)
      ind_balanced = grep('TRUE', ind_balanced)
      child4 = child3[ind_balanced,]
      
      if(nrow(child4) == 0) {
        result[[equation[i]]] = 0
        next    
      }
      
      # If reaction is transport reaction, just pick the same instance substrand and product
      if(class_substrate == class_product) {
        ind_same_sp = which(child4$chebi.x == child4$chebi.y)
        instance_equation_final = instance_equation[ind_balanced][ind_same_sp]
        result[[equation[i]]] = instance_equation_final 
        next
      }
      
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
          message(sprintf('WARNING: For compound %s, target %s have same similarity, so you have to check manually', k, paste(target[max_target], collapse=' and ')))
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
            message(sprintf('WARNING: %s and %s have similarity %s, you should check manually', tmpChild1, tmpChild2, tmpSimilarity))
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
      ind_x = which(child4[[id_col.x]] %in% unique_table[[id_col.x]])
      ind_y = which(child4[[id_col.y]] %in% unique_table[[id_col.y]])
      ind_max = intersect(ind_x, ind_y)
      
      instance_equation_final = instance_equation[ind_balanced][ind_max]
      
      result[[equation[i]]] = instance_equation_final 
    } else {
      message(sprintf('WARNING: %s is not instantiated by this function', equation[i]))
      result[[equation[i]]] = 0
    }
  }
  message(sprintf('# of instantiated generic reaction: %i', length(number_of_instantiated_generic) - no_pair))
  return(result)
}