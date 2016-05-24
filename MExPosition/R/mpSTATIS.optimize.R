mpSTATIS.optimize <- function(data, num.obs, column.design = NULL, num.groups, optimization.option = 'STATIS') 
{ 
    ## Exception
    if(optimization.option != 'None' && optimization.option != 'STATIS' && optimization.option != 'RV_Matrix' && 
        optimization.option != 'Power_STATIS' && optimization.option != 'ANISOSTATIS_Type1' && optimization.option != 'ANISOSTATIS_Type2'
        && optimization.option != 'PTA' && optimization.option != 'CANOSTATIS')
    {	print(paste('WARNING: Optimization option ', optimization.option, ' not recognized. STATIS was set as default'))
        statis.processed <- mpSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = "STATIS")	
    }

    # No optimization
    if(optimization.option == "None")
    {	print(paste('Optimizing using: ',optimization.option))
    	statis.processed <- mpSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = "None")
    }

    #Optimization Option - STATIS
    if(optimization.option == "STATIS")
    {	print(paste('Optimizing using: ',optimization.option))	
	statis.processed <- mpSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = "STATIS")
    }

    #Optimization Option - RV Matrix
    if(optimization.option == 'RV_Matrix')
    {	print(paste('Optimizing using: ',optimization.option))
	statis.processed <- mpSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = 'RV_Matrix')
    }	

    #Optimization Option - Power STATIS
    if(optimization.option == 'STATIS_Power1')
    { 	print(paste('Optimizing using: ',optimization.option))
	statis.processed <- mpSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = 'STATIS_Power1')
    }

    #Optimization Option - ANISOSTATIS Type 1
    if(optimization.option == 'ANISOSTATIS_Type1')
    {	print(paste('Optimizing using: ',optimization.option))
	statis.processed <- mpANISOSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = 'ANISOSTATIS_Type1')
    }

    #Optimization Option - ANISOSTATIS Type 2
    if(optimization.option == 'ANISOSTATIS_Type2')
    {	print(paste('Optimizing using: ',optimization.option))
	statis.processed <- mpANISOSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = 'ANISOSTATIS_Type2')
    }

    #Optimization Option - PTA
    if(optimization.option == 'PTA')
    {   print(paste('Optimizing using: ',optimization.option))
    statis.processed <- mpSTATIS.core(data, num.obs = num.obs, column.design, num.groups = num.groups, optimization.option = 'PTA')
    }

# Results    
{   res.optimize <- list(S=statis.processed$S, RVMatrix = statis.processed$RVMatrix, C = statis.processed$C, ci = statis.processed$ci, 
                    cj = statis.processed$cj, eigs.vector = statis.processed$eigs.vector, eigs = statis.processed$eigs, 
                    fi = statis.processed$fi, tau = statis.processed$tau, alphaWeights = statis.processed$alphaWeights, 
                    
                    compromise = statis.processed$compromise, compromise.ci = statis.processed$compromise.ci, compromise.cj = statis.processed$compromise.cj,
                    compromise.eigs.vector = statis.processed$compromise.eigs.vector, compromise.eigs = statis.processed$compromise.eigs,
                    compromise.fi = statis.processed$compromise.fi, compromise.tau = statis.processed$compromise.tau,
                    
                    alphaWeights = statis.processed$weights, masses = statis.processed$masses, 
                    table.ci = statis.processed$table.ci, table.cj = statis.processed$table.cj, table.eigs = statis.processed$table.eigs, 
                    table.eigs.vector = statis.processed$table.eigs.vector, table.partial.fi.array = statis.processed$table.partial.fi.array,
                    table.loadings = statis.processed$table.loadings, table.fi = statis.processed$table.fi, table.partial.fi = statis.processed$table.partial.fi,
                    table.tau=statis.processed$table.tau)
}

return(res.optimize)   
print('Processing completed')
}