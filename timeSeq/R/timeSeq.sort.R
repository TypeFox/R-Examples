timeSeq.sort = function(genenames, NPDE.ratio, PDE.ratio, table, count) {
	
  Length = length(genenames)
    
## rank the most significant NPDE genes
  NPDE_order = order(-NPDE.ratio)
  NPDE_list = data.frame(genenames = genenames[NPDE_order], ratio = NPDE.ratio[NPDE_order], count = count[NPDE_order])
  table1 = table[NPDE_order, , ]
	
## rank the most significant PDE genes
  PDE_order = order(-PDE.ratio)
  PDE_list = data.frame(genenames = genenames[PDE_order], ratio = PDE.ratio[PDE_order], count = count[PDE_order])
  table2 = table[PDE_order, , ]
  
  out = list(NPDE_list = NPDE_list, 
             PDE_list = PDE_list,
             table1 = table1, 
             table2 = table2)
  out
}

