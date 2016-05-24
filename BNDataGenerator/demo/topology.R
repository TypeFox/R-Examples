##### Start
set.seed(1234)

packages <- c("bnlearn")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

require(bnlearn)
require(BNDataGenerator)

n = 1000

nodes = 10
temp = make_topology(nodes, topology = "Collapse")
target_arcs_mat = temp$arcs_mat
res = BN_Data_Generator(temp$arcs_mat, temp$Probs, n, temp$nodename)
nodename = res$nodename
data = res$data





# Constraint-based algorithms
bn_gs = gs(data)				# the Grow-Shrink(GS)
bn_gs
# coefficients(bn.fit(bn_gs, data))

bn_iamb = iamb(data)			# the Incremental Association(IAMB)
bn_iamb
# coefficients(bn.fit(bn_iamb, data))


# Score-Based Algorithm
bn_hc = hc(data)				# Hill-Climbing (HC)
# bn_hc
# coefficients(bn.fit(bn_hc, data))
round(score(bn_hc, data, type = "bde"), 2)
round(score(bn_hc, data, type = "k2"), 2)
round(score(bn_hc, data, type = "loglik"), 2)
round(score(bn_hc, data, type = "aic"), 2)
round(score(bn_hc, data, type = "bic"), 2)


bn_tabu = tabu(data)			# Tabu search (TABU)
# bn_tabu
# coefficients(bn.fit(bn_tabu, data))
round(score(bn_tabu, data, type = "bde"), 2)
round(score(bn_tabu, data, type = "k2"), 2)
round(score(bn_tabu, data, type = "loglik"), 2)
round(score(bn_tabu, data, type = "aic"), 2)
round(score(bn_tabu, data, type = "bic"), 2)



# Hybsrid algorithms
bn_mmhc = mmhc(data)		# the Max-Min Hill Climbing (MMHC)
# bn_mmhc
# coefficients(bn.fit(bn_mmhc, data))
round(score(bn_mmhc, data, type = "bde"), 2)
round(score(bn_mmhc, data, type = "k2"), 2)
round(score(bn_mmhc, data, type = "loglik"), 2)
round(score(bn_mmhc, data, type = "aic"), 2)
round(score(bn_mmhc, data, type = "bic"), 2)


bn_rsmax2 = rsmax2(data)	# the more general 2-phase Restricted Maximization (RSMAX2)
# bn_rsmax2
# coefficients(bn.fit(bn_rsmax2, data))
round(score(bn_rsmax2, data, type = "bde"), 2)
round(score(bn_rsmax2, data, type = "k2"), 2)
round(score(bn_rsmax2, data, type = "loglik"), 2)
round(score(bn_rsmax2, data, type = "aic"), 2)
round(score(bn_rsmax2, data, type = "bic"), 2)


bn_hc_arcs_mat = fromto_to_mat(bn_hc$arcs, nodename)
bn_tabu_arcs_mat = fromto_to_mat(bn_tabu$arcs, nodename)
bn_mmhc_arcs_mat = fromto_to_mat(bn_mmhc$arcs, nodename)
bn_rsmax2_arcs_mat = fromto_to_mat(bn_rsmax2$arcs, nodename)


## C		# Correct Edges
## M		# Missed Edges
## WO	# Wrongly Oriented Edges
## WC	# Wrongly Connected Edges
C_M_WO_WC_mat = rbind(	C_M_WO_WC(target_arcs_mat, bn_hc_arcs_mat),				# HC
										C_M_WO_WC(target_arcs_mat, bn_tabu_arcs_mat),			# TABU
										C_M_WO_WC(target_arcs_mat, bn_mmhc_arcs_mat),		# MMHC
										C_M_WO_WC(target_arcs_mat, bn_rsmax2_arcs_mat))		# RSMAX2
dimnames(C_M_WO_WC_mat)[[1]] = c("HC", "TABU", "MMHC", "RSMAX2")
dimnames(C_M_WO_WC_mat)[[2]] = c("C", "M", "WO", "WC")
C_M_WO_WC_mat


			
par(mfrow = c(2, 3))

# Constraint-based algorithms
plot(bn_gs, main = "Grow-Shrink")
plot(bn_iamb, main = "Incremental Association")

# Score-Based Algorithm
plot(bn_hc, main = "Hill-Climbing")		
plot(bn_tabu, main = "Tabu search")		

# Hybrid algorithms
plot(bn_mmhc, main = "Max-Min Hill Climbing")
plot(bn_rsmax2, main = "Restricted Maximization")
