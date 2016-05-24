
## ED evaluate (with very few samples)
output <- evaluate.e.ofv.fim(poped.db,ED_samp_size=10)
output$E_ofv

## API evaluate (with very few samples)
output <- evaluate.e.ofv.fim(poped.db,ED_samp_size=10,ofv_calc_type=4)
output$E_ofv

## ED evaluate using Laplace approximation 
tic()
output <- evaluate.e.ofv.fim(poped.db,use_laplace=TRUE)
toc()
output$E_ofv

\dontrun{

  ## ED expected value with more precision. 
  ## Compare time and value to Laplace approximation.
  ## Run a couple of times to see stochasticity of calculation.
  tic()
  e_ofv_mc <- evaluate.e.ofv.fim(poped.db,ED_samp_size=500)
  toc()
  e_ofv_mc$E_ofv
  
  # If you want to get an E(FIM) from the laplace approximation you have to ask for it
  # and it will take more time.
  output <- evaluate.e.ofv.fim(poped.db,use_laplace=TRUE,laplace.fim=TRUE)
  output$E_fim
  
 

}
