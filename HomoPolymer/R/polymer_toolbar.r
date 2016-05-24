env_polymer<-as.environment(1)

 .onAttach<-function(libname, pkgname){
	polymer_toolbar()
}

polymer_toolbar<-function(){
	# bar definition
	polymerbar<-gtkMenuBar()
	# Data Handling item
	DH_menu<-gtkMenu()
	DH_item<-gtkMenuItemNewWithMnemonic(label="_Homopolymer")
	DH_item$setSubmenu(DH_menu)
	polymerbar$append(DH_item)
	DHprep_item<-gtkMenuItemNewWithMnemonic(label="_Prepare")
	gSignalConnect(DHprep_item,"activate",
				function(item){
					if(!exists('recipe.inp',envir=env_polymer)){
						recipe.inp<-data.frame(rep(-1,6),rep(0,6))
					}else{
						recipe.inp<-get('recipe.inp',envir=env_polymer)
					}
					if(!exists('process.inp',envir=env_polymer)){
						process.inp<-c(rep('0',10),rep('FALSE',10))
					}else{
						process.inp<-get('process.inp',envir=env_polymer)
					}
					if(!exists('flow.inp',envir=env_polymer)){
						flow.inp<-list()
						flow.inp<-lapply(1:8,function(i){flow.inp[[i]]<-c(0,0,0)})
					}else{
						flow.inp<-get('flow.inp',envir=env_polymer)
					}
					DBpp<-get('DBpp',envir=env_polymer)
					DBk<-get('DBk',envir=env_polymer)
					out<-prepare(DBpp,DBk,recipe.inp,process.inp,flow.inp,interact=TRUE)
					assign('recipe.inp',out$recipe.inp,envir=env_polymer)
					assign('process.inp',out$process.inp,envir=env_polymer)
					assign('flow.inp',out$flow.inp,envir=env_polymer)
					assign('out.prepare',out[1:(length(out)-3)],envir=env_polymer)
				})
	DH_menu$append(DHprep_item)
	DHevalu_item<-gtkMenuItemNewWithMnemonic(label="_Evaluate")
	gSignalConnect(DHevalu_item,"activate",
				function(item){
					out.prepare<-get('out.prepare',envir=env_polymer)
					out<-integral(out.prepare)
					print('System integration performed correctly',quote=FALSE)
					assign('out.integral',out,envir=env_polymer)
					PlotBasic(out,out.prepare$pars)
				})					
	DH_menu$append(DHevalu_item)
	DHplot_item<-gtkMenuItemNewWithMnemonic(label="_Plot")
	gSignalConnect(DHplot_item,"activate",
					function(item){
						plotdf<-data.frame(
							symbol=c('time','T','X','M','Vl','HA','AM','Mwm','Mnm','Mw','Mn','pH','roout','R','I','S','P','Z','CTA','CCTA',
							'Na','Cp','Kp','Kpp','Ktb','Kd','BN3','BN4'),
							text=c('Time (min)','Temperature (K)','Conversion (-)','Monomer Conc. (mol/L)','Liquid Volume (L)',
							'Undissociated Monomer (mol/L)','Dissociated Monomer (mol/L)','Accumulated Mw (g/mol)','Accumulated Mn (g/mol)',
							'Actual Mw (g/mol)','Actual Mn (g/mol)','pH','Out flow density (kg/L)','Radical Conc (mol/L)','Initiator Conc. (mol/L)',
							'Solvent Conc. (mol/L)','Polymer Conc. (mol/L)','Inhibitor Conc. (mol/L)','Chain Transfer Conc. (mol/L)',
							'Catalytic CA Conc. (mol/L)','Sodium Cations Conc. (mol/L)','Specific heat Cp (cal/mol k)','Propagation Kp (L/mol min)',
							'Eff. Propag. Kpp(L/mol min)','Termination Kt (L/mol min)','Initiator decomp. Kd (1/min)','BN3','BN4')
						)
						out<-get('out.integral',envir=env_polymer)
						PlotChoice(out,plotdf)
					})
	DH_menu$append(DHplot_item)
	DHreset_item<-gtkMenuItemNewWithMnemonic(label="_Reset")
	gSignalConnect(DHreset_item,"activate",
					function(item){
					if(exists('recipe.inp',envir=env_polymer)){
						recipe.inp<-get('recipe.inp',envir=env_polymer)
						rm(recipe.inp,envir=env_polymer)
					}
					if(exists('process.inp',envir=env_polymer)){
						process.inp<-get('process.inp',envir=env_polymer)
						rm(process.inp,envir=env_polymer)
					}
					if(exists('flow.inp',envir=env_polymer)){
						flow.inp<-get('flow.inp',envir=env_polymer)
						rm(flow.inp,envir=env_polymer)
					}
					if(exists('out.prepare',envir=env_polymer)){
						out.prepare<-get('out.prepare',envir=env_polymer)
						rm(out.prepare,envir=env_polymer)
					}
					if(exists('out.integral',envir=env_polymer))
						out.integral<-get('out.integral',envir=env_polymer)
						rm(out.integral,envir=env_polymer)
					})
	DH_menu$append(DHreset_item)

	# build bar
	polymer_window<-gtkWindow(type='GTK_WINDOW_TOPLEVEL')
	polymer_vbox<-gtkVBox()
	polymer_window$add(polymer_vbox)
	polymer_vbox$packStart(polymerbar,FALSE,FALSE)
	polymer_window$setTitle('Polymerization Menubar')
	polymer_window$SetResizable(FALSE)
	polymer_window$Resize(750,20)
}
