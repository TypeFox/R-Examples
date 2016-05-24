# Automatically created wrapper file

AT.run.CPPSC.method <- function( E.MeV.u,
			particle.no,
			fluence.cm2.or.dose.Gy,
			material.no,
			stopping.power.source.no,
			rdd.model,
			rdd.parameters,
			er.model,
			gamma.model,
			gamma.parameters,
			N2,
			fluence.factor,
			write.output,
			shrink.tails,
			shrink.tails.under,
			adjust.N2,
			lethal.events.mode){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2.or.dose.Gy)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	rdd.parameters <- c(rdd.parameters, rep(0, 4 - length(rdd.parameters)))

	gamma.parameters <- c(gamma.parameters, rep(0, 9 - length(gamma.parameters)))

	relative.efficiency <- numeric(1)
	d.check <- numeric(1)
	S.HCP <- numeric(1)
	S.gamma <- numeric(1)
	mean.number.of.tracks.contrib <- numeric(1)
	start.number.of.tracks.contrib <- numeric(1)
	n.convolutions <- numeric(1)
	lower.Jensen.bound <- numeric(1)
	upper.Jensen.bound <- numeric(1)

	res <- .C("AT_run_CPPSC_method_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2.or.dose.Gy = as.single(fluence.cm2.or.dose.Gy),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			rdd.model = as.integer(rdd.model),
			rdd.parameters = as.single(rdd.parameters),
			er.model = as.integer(er.model),
			gamma.model = as.integer(gamma.model),
			gamma.parameters = as.single(gamma.parameters),
			N2 = as.integer(N2),
			fluence.factor = as.single(fluence.factor),
			write.output = as.integer(write.output),
			shrink.tails = as.integer(shrink.tails),
			shrink.tails.under = as.single(shrink.tails.under),
			adjust.N2 = as.integer(adjust.N2),
			lethal.events.mode = as.integer(lethal.events.mode),
			relative.efficiency = as.single(relative.efficiency),
			d.check = as.single(d.check),
			S.HCP = as.single(S.HCP),
			S.gamma = as.single(S.gamma),
			mean.number.of.tracks.contrib = as.single(mean.number.of.tracks.contrib),
			start.number.of.tracks.contrib = as.single(start.number.of.tracks.contrib),
			n.convolutions = as.integer(n.convolutions),
			lower.Jensen.bound = as.single(lower.Jensen.bound),
			upper.Jensen.bound = as.single(upper.Jensen.bound),PACKAGE="libamtrack")

	 return.list <- list(10)
	 return.list[[1]] <- res$N2
	 return.list[[2]] <- res$relative.efficiency
	 return.list[[3]] <- res$d.check
	 return.list[[4]] <- res$S.HCP
	 return.list[[5]] <- res$S.gamma
	 return.list[[6]] <- res$mean.number.of.tracks.contrib
	 return.list[[7]] <- res$start.number.of.tracks.contrib
	 return.list[[8]] <- res$n.convolutions
	 return.list[[9]] <- res$lower.Jensen.bound
	 return.list[[10]] <- res$upper.Jensen.bound
	 names(return.list) <- c("N2","relative.efficiency","d.check","S.HCP","S.gamma","mean.number.of.tracks.contrib","start.number.of.tracks.contrib","n.convolutions","lower.Jensen.bound","upper.Jensen.bound")
	 return(return.list)
}


AT.run.GSM.method <- function( E.MeV.u,
			particle.no,
			fluence.cm2.or.dose.Gy,
			material.no,
			stopping.power.source.no,
			rdd.model,
			rdd.parameters,
			er.model,
			gamma.model,
			gamma.parameters,
			N.runs,
			write.output,
			nX,
			voxel.size.m,
			lethal.events.mode){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2.or.dose.Gy)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	rdd.parameters <- c(rdd.parameters, rep(0, 4 - length(rdd.parameters)))

	gamma.parameters <- c(gamma.parameters, rep(0, 9 - length(gamma.parameters)))

	relative.efficiency <- numeric(1)
	d.check <- numeric(1)
	S.HCP <- numeric(1)
	S.gamma <- numeric(1)
	n.particles <- numeric(1)
	sd.relative.efficiency <- numeric(1)
	sd.d.check <- numeric(1)
	sd.S.HCP <- numeric(1)
	sd.S.gamma <- numeric(1)
	sd.n.particles <- numeric(1)

	res <- .C("AT_run_GSM_method_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2.or.dose.Gy = as.single(fluence.cm2.or.dose.Gy),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			rdd.model = as.integer(rdd.model),
			rdd.parameters = as.single(rdd.parameters),
			er.model = as.integer(er.model),
			gamma.model = as.integer(gamma.model),
			gamma.parameters = as.single(gamma.parameters),
			N.runs = as.integer(N.runs),
			write.output = as.integer(write.output),
			nX = as.integer(nX),
			voxel.size.m = as.single(voxel.size.m),
			lethal.events.mode = as.integer(lethal.events.mode),
			relative.efficiency = as.single(relative.efficiency),
			d.check = as.single(d.check),
			S.HCP = as.single(S.HCP),
			S.gamma = as.single(S.gamma),
			n.particles = as.single(n.particles),
			sd.relative.efficiency = as.single(sd.relative.efficiency),
			sd.d.check = as.single(sd.d.check),
			sd.S.HCP = as.single(sd.S.HCP),
			sd.S.gamma = as.single(sd.S.gamma),
			sd.n.particles = as.single(sd.n.particles),PACKAGE="libamtrack")

	 return.list <- list(10)
	 return.list[[1]] <- res$relative.efficiency
	 return.list[[2]] <- res$d.check
	 return.list[[3]] <- res$S.HCP
	 return.list[[4]] <- res$S.gamma
	 return.list[[5]] <- res$n.particles
	 return.list[[6]] <- res$sd.relative.efficiency
	 return.list[[7]] <- res$sd.d.check
	 return.list[[8]] <- res$sd.S.HCP
	 return.list[[9]] <- res$sd.S.gamma
	 return.list[[10]] <- res$sd.n.particles
	 names(return.list) <- c("relative.efficiency","d.check","S.HCP","S.gamma","n.particles","sd.relative.efficiency","sd.d.check","sd.S.HCP","sd.S.gamma","sd.n.particles")
	 return(return.list)
}


AT.run.IGK.method <- function( E.MeV.u,
			particle.no,
			fluence.cm2.or.dose.Gy,
			material.no,
			stopping.power.source.no,
			rdd.model,
			rdd.parameters,
			er.model,
			gamma.model,
			gamma.parameters,
			saturation.cross.section.factor,
			write.output){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2.or.dose.Gy)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	rdd.parameters <- c(rdd.parameters, rep(0, 4 - length(rdd.parameters)))

	gamma.parameters <- c(gamma.parameters, rep(0, 9 - length(gamma.parameters)))

	relative.efficiency <- numeric(1)
	S.HCP <- numeric(1)
	S.gamma <- numeric(1)
	sI.cm2 <- numeric(1)
	gamma.dose.Gy <- numeric(1)
	P.I <- numeric(1)
	P.g <- numeric(1)

	res <- .C("AT_run_IGK_method_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2.or.dose.Gy = as.single(fluence.cm2.or.dose.Gy),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			rdd.model = as.integer(rdd.model),
			rdd.parameters = as.single(rdd.parameters),
			er.model = as.integer(er.model),
			gamma.model = as.integer(gamma.model),
			gamma.parameters = as.single(gamma.parameters),
			saturation.cross.section.factor = as.single(saturation.cross.section.factor),
			write.output = as.integer(write.output),
			relative.efficiency = as.single(relative.efficiency),
			S.HCP = as.single(S.HCP),
			S.gamma = as.single(S.gamma),
			sI.cm2 = as.single(sI.cm2),
			gamma.dose.Gy = as.single(gamma.dose.Gy),
			P.I = as.single(P.I),
			P.g = as.single(P.g),PACKAGE="libamtrack")

	 return.list <- list(7)
	 return.list[[1]] <- res$relative.efficiency
	 return.list[[2]] <- res$S.HCP
	 return.list[[3]] <- res$S.gamma
	 return.list[[4]] <- res$sI.cm2
	 return.list[[5]] <- res$gamma.dose.Gy
	 return.list[[6]] <- res$P.I
	 return.list[[7]] <- res$P.g
	 names(return.list) <- c("relative.efficiency","S.HCP","S.gamma","sI.cm2","gamma.dose.Gy","P.I","P.g")
	 return(return.list)
}


AT.set.user.material.from.composition <- function( density.g.cm3,
			Z,
			A,
			weight.fraction){

	n	<- length(A)
	if(n != length(Z)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(weight.fraction)){cat("Array size mismatch for 'n'!\n")
		return}

	status <- numeric(1)

	res <- .C("AT_set_user_material_from_composition_R", 
			n = as.integer(n),
			density.g.cm3 = as.single(density.g.cm3),
			A = as.integer(A),
			Z = as.integer(Z),
			weight.fraction = as.single(weight.fraction),
			status = as.integer(status),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$status
	 names(return.list) <- c("status")
	 return(return.list)
}


AT.set.user.material <- function( density.g.cm3,
			I.eV,
			average.A,
			average.Z){

	status <- numeric(1)

	res <- .C("AT_set_user_material_R", 
			density.g.cm3 = as.single(density.g.cm3),
			I.eV = as.single(I.eV),
			average.A = as.single(average.A),
			average.Z = as.single(average.Z),
			status = as.integer(status),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$status
	 names(return.list) <- c("status")
	 return(return.list)
}


AT.I.eV.from.composition <- function( Z,
			A,
			weight.fraction){

	n	<- length(Z)
	if(n != length(A)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(weight.fraction)){cat("Array size mismatch for 'n'!\n")
		return}

	I.eV <- numeric(1)

	res <- .C("AT_I_eV_from_composition_R", 
			n = as.integer(n),
			Z = as.integer(Z),
			A = as.integer(A),
			weight.fraction = as.single(weight.fraction),
			I.eV = as.single(I.eV),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$I.eV
	 names(return.list) <- c("I.eV")
	 return(return.list)
}


AT.effective.Z.from.composition <- function( Z,
			weight.fraction,
			electron.densities.cm3,
			exponent){

	n	<- length(Z)
	if(n != length(weight.fraction)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(electron.densities.cm3)){cat("Array size mismatch for 'n'!\n")
		return}

	effective.Z <- numeric(1)

	res <- .C("AT_effective_Z_from_composition_R", 
			n = as.integer(n),
			Z = as.integer(Z),
			weight.fraction = as.single(weight.fraction),
			electron.densities.cm3 = as.single(electron.densities.cm3),
			exponent = as.single(exponent),
			effective.Z = as.single(effective.Z),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$effective.Z
	 names(return.list) <- c("effective.Z")
	 return(return.list)
}


AT.average.Z.from.composition <- function( Z,
			weight.fraction){

	n	<- length(Z)
	if(n != length(weight.fraction)){cat("Array size mismatch for 'n'!\n")
		return}

	average.Z <- numeric(1)

	res <- .C("AT_average_Z_from_composition_R", 
			n = as.integer(n),
			Z = as.integer(Z),
			weight.fraction = as.single(weight.fraction),
			average.Z = as.single(average.Z),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$average.Z
	 names(return.list) <- c("average.Z")
	 return(return.list)
}


AT.average.A.from.composition <- function( A,
			weight.fraction){

	n	<- length(A)
	if(n != length(weight.fraction)){cat("Array size mismatch for 'n'!\n")
		return}

	average.A <- numeric(1)

	res <- .C("AT_average_A_from_composition_R", 
			n = as.integer(n),
			A = as.integer(A),
			weight.fraction = as.single(weight.fraction),
			average.A = as.single(average.A),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$average.A
	 names(return.list) <- c("average.A")
	 return(return.list)
}


AT.electron.density.m3.from.composition <- function( density.g.cm3,
			Z,
			A,
			weight.fraction){

	n	<- length(Z)
	if(n != length(A)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(weight.fraction)){cat("Array size mismatch for 'n'!\n")
		return}

	electron.density.m3 <- numeric(1)

	res <- .C("AT_electron_density_m3_from_composition_R", 
			n = as.integer(n),
			density.g.cm3 = as.single(density.g.cm3),
			Z = as.integer(Z),
			A = as.integer(A),
			weight.fraction = as.single(weight.fraction),
			electron.density.m3 = as.single(electron.density.m3),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$electron.density.m3
	 names(return.list) <- c("electron.density.m3")
	 return(return.list)
}


AT.electron.density.m3 <- function( density.g.cm3,
			average.A,
			average.Z){

	n	<- length(density.g.cm3)
	if(n != length(average.Z)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(average.A)){cat("Array size mismatch for 'n'!\n")
		return}

	electron.density.m3 <- numeric(n)

	res <- .C("AT_electron_density_m3_multi_R", 
			n = as.integer(n),
			density.g.cm3 = as.single(density.g.cm3),
			average.Z = as.single(average.Z),
			average.A = as.single(average.A),
			electron.density.m3 = as.single(electron.density.m3),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$electron.density.m3
	 names(return.list) <- c("electron.density.m3")
	 return(return.list)
}


AT.electron.density.m3.from.material.no <- function( material.no){

	n	<- length(material.no)
	electron.density.m3 <- numeric(n)

	res <- .C("AT_electron_density_m3_from_material_no_multi_R", 
			n = as.integer(n),
			material.no = as.integer(material.no),
			electron.density.m3 = as.single(electron.density.m3),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$electron.density.m3
	 names(return.list) <- c("electron.density.m3")
	 return(return.list)
}


AT.get.materials.data <- function( material.no){

	number.of.materials	<- length(material.no)
	density.g.cm3 <- numeric(number.of.materials)
	I.eV <- numeric(number.of.materials)
	alpha.g.cm2.MeV <- numeric(number.of.materials)
	p.MeV <- numeric(number.of.materials)
	m.g.cm2 <- numeric(number.of.materials)
	average.A <- numeric(number.of.materials)
	average.Z <- numeric(number.of.materials)

	res <- .C("AT_get_materials_data_R", 
			number.of.materials = as.integer(number.of.materials),
			material.no = as.integer(material.no),
			density.g.cm3 = as.single(density.g.cm3),
			I.eV = as.single(I.eV),
			alpha.g.cm2.MeV = as.single(alpha.g.cm2.MeV),
			p.MeV = as.single(p.MeV),
			m.g.cm2 = as.single(m.g.cm2),
			average.A = as.single(average.A),
			average.Z = as.single(average.Z),PACKAGE="libamtrack")

	 return.list <- list(7)
	 return.list[[1]] <- res$density.g.cm3
	 return.list[[2]] <- res$I.eV
	 return.list[[3]] <- res$alpha.g.cm2.MeV
	 return.list[[4]] <- res$p.MeV
	 return.list[[5]] <- res$m.g.cm2
	 return.list[[6]] <- res$average.A
	 return.list[[7]] <- res$average.Z
	 names(return.list) <- c("density.g.cm3","I.eV","alpha.g.cm2.MeV","p.MeV","m.g.cm2","average.A","average.Z")
	 return(return.list)
}


AT.nuclear.spin.from.particle.no <- function( particle.no){

	n	<- length(particle.no)
	I <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_nuclear_spin_from_particle_no_multi_R", 
			n = as.integer(n),
			particle.no = as.integer(particle.no),
			I = as.single(I),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$I
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("I", "returnValue")
	 return(return.list)
}


AT.Z.from.particle.no <- function( particle.no){

	n	<- length(particle.no)
	Z <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_Z_from_particle_no_R", 
			n = as.integer(n),
			particle.no = as.integer(particle.no),
			Z = as.integer(Z),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$Z
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("Z", "returnValue")
	 return(return.list)
}


AT.atomic.weight.from.Z <- function( Z){

	n	<- length(Z)
	atomic.weight <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_atomic_weight_from_Z_R", 
			n = as.integer(n),
			Z = as.integer(Z),
			atomic.weight = as.single(atomic.weight),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$atomic.weight
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("atomic.weight", "returnValue")
	 return(return.list)
}


AT.A.from.particle.no <- function( particle.no){

	n	<- length(particle.no)
	A <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_A_from_particle_no_R", 
			n = as.integer(n),
			particle.no = as.integer(particle.no),
			A = as.integer(A),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$A
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("A", "returnValue")
	 return(return.list)
}


AT.particle.no.from.Z.and.A <- function( Z,
			A){

	n	<- length(Z)
	if(n != length(A)){cat("Array size mismatch for 'n'!\n")
		return}

	particle.no <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_particle_no_from_Z_and_A_R", 
			n = as.integer(n),
			Z = as.integer(Z),
			A = as.integer(A),
			particle.no = as.integer(particle.no),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$particle.no
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("particle.no", "returnValue")
	 return(return.list)
}


AT.CSDA.range.g.cm2 <- function( E.initial.MeV.u,
			E.final.MeV.u,
			particle.no,
			material.no){

	n	<- length(E.initial.MeV.u)
	if(n != length(E.final.MeV.u)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	CSDA.range.cm2.g <- numeric(n)

	res <- .C("AT_CSDA_range_g_cm2_multi_R", 
			n = as.integer(n),
			E.initial.MeV.u = as.single(E.initial.MeV.u),
			E.final.MeV.u = as.single(E.final.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			CSDA.range.cm2.g = as.single(CSDA.range.cm2.g),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$CSDA.range.cm2.g
	 names(return.list) <- c("CSDA.range.cm2.g")
	 return(return.list)
}


AT.max.electron.ranges.m <- function( E.MeV.u,
			material.no,
			er.model){

	number.of.particles	<- length(E.MeV.u)
	max.electron.range.m <- numeric(number.of.particles)

	res <- .C("AT_max_electron_ranges_m_R", 
			number.of.particles = as.integer(number.of.particles),
			E.MeV.u = as.single(E.MeV.u),
			material.no = as.integer(material.no),
			er.model = as.integer(er.model),
			max.electron.range.m = as.single(max.electron.range.m),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$max.electron.range.m
	 names(return.list) <- c("max.electron.range.m")
	 return(return.list)
}


AT.mean.energy.loss.keV <- function( E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	returnValue = numeric(1)

	res <- .C("AT_mean_energy_loss_keV_R", 
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.xi.keV <- function( E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	returnValue = numeric(1)

	res <- .C("AT_xi_keV_R", 
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.kappa <- function( E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(slab.thickness.um)){cat("Array size mismatch for 'n'!\n")
		return}

	kappa <- numeric(n)

	res <- .C("AT_kappa_multi_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			kappa = as.single(kappa),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$kappa
	 names(return.list) <- c("kappa")
	 return(return.list)
}


AT.Landau.PDF <- function( lambda.landau){

	n	<- length(lambda.landau)
	density <- numeric(n)

	res <- .C("AT_Landau_PDF_R", 
			n = as.integer(n),
			lambda.landau = as.single(lambda.landau),
			density = as.single(density),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$density
	 names(return.list) <- c("density")
	 return(return.list)
}


AT.Landau.IDF <- function( rnd){

	n	<- length(rnd)
	lambda.landau <- numeric(n)

	res <- .C("AT_Landau_IDF_R", 
			n = as.integer(n),
			rnd = as.single(rnd),
			lambda.landau = as.single(lambda.landau),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$lambda.landau
	 names(return.list) <- c("lambda.landau")
	 return(return.list)
}


AT.lambda.landau.from.energy.loss <- function( energy.loss.keV,
			E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(energy.loss.keV)
	lambda.landau <- numeric(n)

	res <- .C("AT_lambda_landau_from_energy_loss_multi_R", 
			n = as.integer(n),
			energy.loss.keV = as.single(energy.loss.keV),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			lambda.landau = as.single(lambda.landau),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$lambda.landau
	 names(return.list) <- c("lambda.landau")
	 return(return.list)
}


AT.lambda.mean <- function( E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(slab.thickness.um)){cat("Array size mismatch for 'n'!\n")
		return}

	lambda.mean <- numeric(n)

	res <- .C("AT_lambda_mean_multi_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			lambda.mean = as.single(lambda.mean),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$lambda.mean
	 names(return.list) <- c("lambda.mean")
	 return(return.list)
}


AT.lambda.max <- function( E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(slab.thickness.um)){cat("Array size mismatch for 'n'!\n")
		return}

	lambda.max <- numeric(n)

	res <- .C("AT_lambda_max_multi_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			lambda.max = as.single(lambda.max),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$lambda.max
	 names(return.list) <- c("lambda.max")
	 return(return.list)
}


AT.energy.loss.from.lambda.landau <- function( lambda.landau,
			E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(lambda.landau)
	if(n != length(E.MeV.u)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(slab.thickness.um)){cat("Array size mismatch for 'n'!\n")
		return}

	energy.loss.keV <- numeric(n)

	res <- .C("AT_energy_loss_from_lambda_landau_multi_R", 
			n = as.integer(n),
			lambda.landau = as.single(lambda.landau),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			energy.loss.keV = as.single(energy.loss.keV),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$energy.loss.keV
	 names(return.list) <- c("energy.loss.keV")
	 return(return.list)
}


AT.Vavilov.PDF <- function( lambda.vavilov,
			kappa,
			beta){

	n	<- length(lambda.vavilov)
	density <- numeric(n)

	res <- .C("AT_Vavilov_PDF_R", 
			n = as.integer(n),
			lambda.vavilov = as.single(lambda.vavilov),
			kappa = as.single(kappa),
			beta = as.single(beta),
			density = as.single(density),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$density
	 names(return.list) <- c("density")
	 return(return.list)
}


AT.Vavilov.IDF <- function( rnd,
			kappa,
			beta){

	n	<- length(rnd)
	if(n != length(kappa)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(beta)){cat("Array size mismatch for 'n'!\n")
		return}

	lambda.vavilov <- numeric(n)

	res <- .C("AT_Vavilov_IDF_R", 
			n = as.integer(n),
			rnd = as.single(rnd),
			kappa = as.single(kappa),
			beta = as.single(beta),
			lambda.vavilov = as.single(lambda.vavilov),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$lambda.vavilov
	 names(return.list) <- c("lambda.vavilov")
	 return(return.list)
}


AT.lambda.vavilov.from.energy.loss <- function( energy.loss.keV,
			E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(energy.loss.keV)
	lambda.vavilov <- numeric(n)

	res <- .C("AT_lambda_vavilov_from_energy_loss_multi_R", 
			n = as.integer(n),
			energy.loss.keV = as.single(energy.loss.keV),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			lambda.vavilov = as.single(lambda.vavilov),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$lambda.vavilov
	 names(return.list) <- c("lambda.vavilov")
	 return(return.list)
}


AT.energy.loss.from.lambda.vavilov <- function( lambda.vavilov,
			E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(lambda.vavilov)
	if(n != length(E.MeV.u)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(slab.thickness.um)){cat("Array size mismatch for 'n'!\n")
		return}

	energy.loss.keV <- numeric(n)

	res <- .C("AT_energy_loss_from_lambda_vavilov_multi_R", 
			n = as.integer(n),
			lambda.vavilov = as.single(lambda.vavilov),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			energy.loss.keV = as.single(energy.loss.keV),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$energy.loss.keV
	 names(return.list) <- c("energy.loss.keV")
	 return(return.list)
}


AT.Gauss.PDF <- function( lambda.gauss){

	n	<- length(lambda.gauss)
	density <- numeric(n)

	res <- .C("AT_Gauss_PDF_R", 
			n = as.integer(n),
			lambda.gauss = as.single(lambda.gauss),
			density = as.single(density),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$density
	 names(return.list) <- c("density")
	 return(return.list)
}


AT.Gauss.IDF <- function( rnd){

	n	<- length(rnd)
	lambda.gauss <- numeric(n)

	res <- .C("AT_Gauss_IDF_R", 
			n = as.integer(n),
			rnd = as.single(rnd),
			lambda.gauss = as.single(lambda.gauss),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$lambda.gauss
	 names(return.list) <- c("lambda.gauss")
	 return(return.list)
}


AT.energy.loss.from.lambda.gauss <- function( lambda.gauss,
			E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.um){

	n	<- length(lambda.gauss)
	if(n != length(E.MeV.u)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(slab.thickness.um)){cat("Array size mismatch for 'n'!\n")
		return}

	energy.loss.keV <- numeric(n)

	res <- .C("AT_energy_loss_from_lambda_gauss_multi_R", 
			n = as.integer(n),
			lambda.gauss = as.single(lambda.gauss),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.um = as.single(slab.thickness.um),
			energy.loss.keV = as.single(energy.loss.keV),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$energy.loss.keV
	 names(return.list) <- c("energy.loss.keV")
	 return(return.list)
}


AT.beta.from.E <- function( E.MeV.u){

	n	<- length(E.MeV.u)
	beta <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_beta_from_E_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			beta = as.single(beta),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$beta
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("beta", "returnValue")
	 return(return.list)
}


AT.E.from.beta <- function( beta){

	n	<- length(beta)
	E.MeV.u <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_E_from_beta_R", 
			n = as.integer(n),
			beta = as.single(beta),
			E.MeV.u = as.single(E.MeV.u),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$E.MeV.u
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("E.MeV.u", "returnValue")
	 return(return.list)
}


AT.E.MeV.u.from.momentum.MeV.c.u <- function( momentum.MeV.c.u){

	n	<- length(momentum.MeV.c.u)
	E.MeV.u <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_E_MeV_u_from_momentum_MeV_c_u_R", 
			n = as.integer(n),
			momentum.MeV.c.u = as.single(momentum.MeV.c.u),
			E.MeV.u = as.single(E.MeV.u),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$E.MeV.u
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("E.MeV.u", "returnValue")
	 return(return.list)
}


AT.gamma.from.E <- function( E.MeV.u){

	n	<- length(E.MeV.u)
	gamma <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_gamma_from_E_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			gamma = as.single(gamma),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$gamma
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("gamma", "returnValue")
	 return(return.list)
}


AT.energy.straggling.MeV2.cm2.g <- function( E.MeV.u,
			particle.no,
			material.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	dsE2dz.MeV2.cm2.g <- numeric(n)

	res <- .C("AT_energy_straggling_MeV2_cm2_g_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			dsE2dz.MeV2.cm2.g = as.single(dsE2dz.MeV2.cm2.g),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$dsE2dz.MeV2.cm2.g
	 names(return.list) <- c("dsE2dz.MeV2.cm2.g")
	 return(return.list)
}


AT.energy.straggling.after.slab.E.MeV.u <- function( E.MeV.u,
			particle.no,
			material.no,
			slab.thickness.m,
			initial.sigma.E.MeV.u){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(initial.sigma.E.MeV.u)){cat("Array size mismatch for 'n'!\n")
		return}

	sigma.E.MeV.u <- numeric(n)

	res <- .C("AT_energy_straggling_after_slab_E_MeV_u_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			slab.thickness.m = as.single(slab.thickness.m),
			initial.sigma.E.MeV.u = as.single(initial.sigma.E.MeV.u),
			sigma.E.MeV.u = as.single(sigma.E.MeV.u),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$sigma.E.MeV.u
	 names(return.list) <- c("sigma.E.MeV.u")
	 return(return.list)
}


AT.effective.charge.from.E.MeV.u <- function( E.MeV.u,
			particle.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	effective.charge <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_effective_charge_from_E_MeV_u_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			effective.charge = as.single(effective.charge),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$effective.charge
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("effective.charge", "returnValue")
	 return(return.list)
}


AT.max.E.transfer.MeV.new <- function( E.MeV.u,
			A){

	n	<- length(E.MeV.u)
	if(n != length(A)){cat("Array size mismatch for 'n'!\n")
		return}

	max.E.transfer.MeV <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_max_E_transfer_MeV_new_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			A = as.integer(A),
			max.E.transfer.MeV = as.single(max.E.transfer.MeV),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$max.E.transfer.MeV
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("max.E.transfer.MeV", "returnValue")
	 return(return.list)
}


AT.max.E.transfer.MeV <- function( E.MeV.u){

	n	<- length(E.MeV.u)
	max.E.transfer.MeV <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_max_E_transfer_MeV_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			max.E.transfer.MeV = as.single(max.E.transfer.MeV),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$max.E.transfer.MeV
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("max.E.transfer.MeV", "returnValue")
	 return(return.list)
}


AT.momentum.MeV.c.u.from.E.MeV.u <- function( E.MeV.u){

	n	<- length(E.MeV.u)
	momentum.MeV.c <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_momentum_MeV_c_u_from_E_MeV_u_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			momentum.MeV.c = as.single(momentum.MeV.c),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$momentum.MeV.c
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("momentum.MeV.c", "returnValue")
	 return(return.list)
}


AT.dose.Gy.from.fluence.cm2 <- function( E.MeV.u,
			particle.no,
			fluence.cm2,
			material.no,
			stopping.power.source.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(fluence.cm2)){cat("Array size mismatch for 'n'!\n")
		return}

	dose.Gy <- numeric(n)

	res <- .C("AT_dose_Gy_from_fluence_cm2_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2 = as.single(fluence.cm2),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			dose.Gy = as.single(dose.Gy),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$dose.Gy
	 names(return.list) <- c("dose.Gy")
	 return(return.list)
}


AT.fluence.cm2.from.dose.Gy <- function( E.MeV.u,
			particle.no,
			D.Gy,
			material.no,
			stopping.power.source.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	if(n != length(D.Gy)){cat("Array size mismatch for 'n'!\n")
		return}

	fluence.cm2 <- numeric(n)

	res <- .C("AT_fluence_cm2_from_dose_Gy_R", 
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			D.Gy = as.single(D.Gy),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			fluence.cm2 = as.single(fluence.cm2),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$fluence.cm2
	 names(return.list) <- c("fluence.cm2")
	 return(return.list)
}


AT.beam.par.physical.to.technical <- function( fluence.cm2,
			sigma.cm){

	n	<- length(fluence.cm2)
	if(n != length(sigma.cm)){cat("Array size mismatch for 'n'!\n")
		return}

	N <- numeric(n)
	FWHM.mm <- numeric(n)

	res <- .C("AT_beam_par_physical_to_technical_R", 
			n = as.integer(n),
			fluence.cm2 = as.single(fluence.cm2),
			sigma.cm = as.single(sigma.cm),
			N = as.single(N),
			FWHM.mm = as.single(FWHM.mm),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$N
	 return.list[[2]] <- res$FWHM.mm
	 names(return.list) <- c("N","FWHM.mm")
	 return(return.list)
}


AT.beam.par.technical.to.physical <- function( N,
			FWHM.mm){

	n	<- length(N)
	if(n != length(FWHM.mm)){cat("Array size mismatch for 'n'!\n")
		return}

	fluence.cm2 <- numeric(n)
	sigma.cm <- numeric(n)

	res <- .C("AT_beam_par_technical_to_physical_R", 
			n = as.integer(n),
			N = as.single(N),
			FWHM.mm = as.single(FWHM.mm),
			fluence.cm2 = as.single(fluence.cm2),
			sigma.cm = as.single(sigma.cm),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$fluence.cm2
	 return.list[[2]] <- res$sigma.cm
	 names(return.list) <- c("fluence.cm2","sigma.cm")
	 return(return.list)
}


AT.total.D.Gy <- function( E.MeV.u,
			particle.no,
			fluence.cm2,
			material.no,
			stopping.power.source.no){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_total_D_Gy_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2 = as.single(fluence.cm2),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.total.fluence.cm2 <- function( E.MeV.u,
			particle.no,
			D.Gy,
			material.no,
			stopping.power.source.no){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(D.Gy)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_total_fluence_cm2_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			D.Gy = as.single(D.Gy),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.fluence.weighted.E.MeV.u <- function( E.MeV.u,
			fluence.cm2){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(fluence.cm2)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_fluence_weighted_E_MeV_u_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			fluence.cm2 = as.single(fluence.cm2),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.dose.weighted.E.MeV.u <- function( E.MeV.u,
			particle.no,
			fluence.cm2,
			material.no,
			stopping.power.source.no){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_dose_weighted_E_MeV_u_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2 = as.single(fluence.cm2),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.fluence.weighted.LET.MeV.cm2.g <- function( E.MeV.u,
			particle.no,
			fluence.cm2,
			material.no,
			stopping.power.source.no){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_fluence_weighted_LET_MeV_cm2_g_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2 = as.single(fluence.cm2),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.dose.weighted.LET.MeV.cm2.g <- function( E.MeV.u,
			particle.no,
			fluence.cm2,
			material.no,
			stopping.power.source.no){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_dose_weighted_LET_MeV_cm2_g_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2 = as.single(fluence.cm2),
			material.no = as.integer(material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.stopping.power.ratio <- function( E.MeV.u,
			particle.no,
			fluence.cm2,
			material.no,
			reference.material.no,
			stopping.power.source.no){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_stopping_power_ratio_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2 = as.single(fluence.cm2),
			material.no = as.integer(material.no),
			reference.material.no = as.integer(reference.material.no),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.mean.number.of.tracks.contrib <- function( E.MeV.u,
			particle.no,
			fluence.cm2,
			material.no,
			er.model,
			stopping.power.source.no){

	number.of.field.components	<- length(E.MeV.u)
	if(number.of.field.components != length(particle.no)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	if(number.of.field.components != length(fluence.cm2)){cat("Array size mismatch for 'number.of.field.components'!\n")
		return}

	returnValue = numeric(1)

	res <- .C("AT_mean_number_of_tracks_contrib_R", 
			number.of.field.components = as.integer(number.of.field.components),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			fluence.cm2 = as.single(fluence.cm2),
			material.no = as.integer(material.no),
			er.model = as.integer(er.model),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			returnValue = as.single(returnValue),PACKAGE="libamtrack")

	 return.list <- list(1)
	 return.list[[1]] <- res$returnValue
	 names(return.list) <- c("returnValue")
	 return(return.list)
}


AT.Rutherford.SDCS <- function( E.MeV.u,
			particle.no,
			material.no,
			T.MeV){

	n	<- length(T.MeV)
	dsdT.m2.MeV <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_Rutherford_SDCS_R", 
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			n = as.integer(n),
			T.MeV = as.single(T.MeV),
			dsdT.m2.MeV = as.single(dsdT.m2.MeV),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$dsdT.m2.MeV
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("dsdT.m2.MeV", "returnValue")
	 return(return.list)
}


AT.r.RDD.m <- function( D.RDD.Gy,
			E.MeV.u,
			particle.no,
			material.no,
			rdd.model,
			rdd.parameter,
			er.model,
			stopping.power.source.no){

	n	<- length(D.RDD.Gy)
	rdd.parameter <- c(rdd.parameter, rep(0, 4 - length(rdd.parameter)))

	r.RDD.m <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_r_RDD_m_R", 
			n = as.integer(n),
			D.RDD.Gy = as.single(D.RDD.Gy),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			rdd.model = as.integer(rdd.model),
			rdd.parameter = as.single(rdd.parameter),
			er.model = as.integer(er.model),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			r.RDD.m = as.single(r.RDD.m),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$r.RDD.m
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("r.RDD.m", "returnValue")
	 return(return.list)
}


AT.D.RDD.Gy <- function( r.m,
			E.MeV.u,
			particle.no,
			material.no,
			rdd.model,
			rdd.parameter,
			er.model,
			stopping.power.source.no){

	n	<- length(r.m)
	rdd.parameter <- c(rdd.parameter, rep(0, 4 - length(rdd.parameter)))

	D.RDD.Gy <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_D_RDD_Gy_R", 
			n = as.integer(n),
			r.m = as.single(r.m),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			rdd.model = as.integer(rdd.model),
			rdd.parameter = as.single(rdd.parameter),
			er.model = as.integer(er.model),
			stopping.power.source.no = as.integer(stopping.power.source.no),
			D.RDD.Gy = as.single(D.RDD.Gy),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$D.RDD.Gy
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("D.RDD.Gy", "returnValue")
	 return(return.list)
}


AT.Mass.Stopping.Power <- function( stopping.power.source,
			E.MeV.u,
			particle.no,
			material.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	stopping.power.MeV.cm2.g <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_Mass_Stopping_Power_R", 
			stopping.power.source = as.character(stopping.power.source),
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			stopping.power.MeV.cm2.g = as.single(stopping.power.MeV.cm2.g),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$stopping.power.MeV.cm2.g
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("stopping.power.MeV.cm2.g", "returnValue")
	 return(return.list)
}


AT.Stopping.Power <- function( stopping.power.source,
			E.MeV.u,
			particle.no,
			material.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	stopping.power.keV.um <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_Stopping_Power_R", 
			stopping.power.source = as.character(stopping.power.source),
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			stopping.power.keV.um = as.single(stopping.power.keV.um),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$stopping.power.keV.um
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("stopping.power.keV.um", "returnValue")
	 return(return.list)
}


AT.Mass.Stopping.Power.with.no <- function( stopping.power.source.no,
			E.MeV.u,
			particle.no,
			material.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	stopping.power.MeV.cm2.g <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_Mass_Stopping_Power_with_no_R", 
			stopping.power.source.no = as.integer(stopping.power.source.no),
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			stopping.power.MeV.cm2.g = as.single(stopping.power.MeV.cm2.g),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$stopping.power.MeV.cm2.g
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("stopping.power.MeV.cm2.g", "returnValue")
	 return(return.list)
}


AT.Stopping.Power.with.no <- function( stopping.power.source.no,
			E.MeV.u,
			particle.no,
			material.no){

	n	<- length(E.MeV.u)
	if(n != length(particle.no)){cat("Array size mismatch for 'n'!\n")
		return}

	stopping.power.keV.um <- numeric(n)
	returnValue = numeric(1)

	res <- .C("AT_Stopping_Power_with_no_R", 
			stopping.power.source.no = as.integer(stopping.power.source.no),
			n = as.integer(n),
			E.MeV.u = as.single(E.MeV.u),
			particle.no = as.integer(particle.no),
			material.no = as.integer(material.no),
			stopping.power.keV.um = as.single(stopping.power.keV.um),
			returnValue = as.integer(returnValue),PACKAGE="libamtrack")

	 return.list <- list(2)
	 return.list[[1]] <- res$stopping.power.keV.um
	 return.list[[2]] <- res$returnValue
	 names(return.list) <- c("stopping.power.keV.um", "returnValue")
	 return(return.list)
}


AT.translate.dose.into.DSB.distribution <- function( f.d.Gy,
			f.dd.Gy,
			f,
			enhancement.factor,
			DSB.per.Gy.per.domain,
			domains.per.nucleus,
			write.output){

	n.bins.f	<- length(f.d.Gy)
	if(n.bins.f != length(f.dd.Gy)){cat("Array size mismatch for 'n.bins.f'!\n")
		return}

	if(n.bins.f != length(f)){cat("Array size mismatch for 'n.bins.f'!\n")
		return}

	if(n.bins.f != length(enhancement.factor)){cat("Array size mismatch for 'n.bins.f'!\n")
		return}

	total.pDSBs <- numeric(1)
	total.nDSBs <- numeric(1)
	number.of.iDSBs <- numeric(1)
	number.of.cDSBs <- numeric(1)
	avg.number.of.DSBs.in.cDSBs <- numeric(1)

	res <- .C("AT_translate_dose_into_DSB_distribution_R", 
			n.bins.f = as.integer(n.bins.f),
			f.d.Gy = as.single(f.d.Gy),
			f.dd.Gy = as.single(f.dd.Gy),
			f = as.single(f),
			enhancement.factor = as.single(enhancement.factor),
			DSB.per.Gy.per.domain = as.single(DSB.per.Gy.per.domain),
			domains.per.nucleus = as.integer(domains.per.nucleus),
			write.output = as.integer(write.output),
			total.pDSBs = as.single(total.pDSBs),
			total.nDSBs = as.single(total.nDSBs),
			number.of.iDSBs = as.single(number.of.iDSBs),
			number.of.cDSBs = as.single(number.of.cDSBs),
			avg.number.of.DSBs.in.cDSBs = as.single(avg.number.of.DSBs.in.cDSBs),PACKAGE="libamtrack")

	 return.list <- list(5)
	 return.list[[1]] <- res$total.pDSBs
	 return.list[[2]] <- res$total.nDSBs
	 return.list[[3]] <- res$number.of.iDSBs
	 return.list[[4]] <- res$number.of.cDSBs
	 return.list[[5]] <- res$avg.number.of.DSBs.in.cDSBs
	 names(return.list) <- c("total.pDSBs","total.nDSBs","number.of.iDSBs","number.of.cDSBs","avg.number.of.DSBs.in.cDSBs")
	 return(return.list)
}


