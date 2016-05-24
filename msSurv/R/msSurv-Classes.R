####################################################################
## Methods/Class Definitions
####################################################################


setClassUnion("AorN", c("array", "NULL"))

setClass("msSurv",representation(tree="graphNEL",
                                 ns="numeric",
                                 et="numeric",
                                 pos.trans="character",
                                 nt.states="character",
                                 dNs="array",
                                 Ys="array",
                                 sum_dNs="array",
                                 dNs.K="array",
                                 Ys.K="array",
                                 sum_dNs.K="array",
                                 ps="array",
                                 AJs="array",
				 I.dA="array",
                                 cov.AJs="AorN",
                                 var.sop="AorN",
				 cov.dA="AorN",
                                 Fnorm="AorN",
                                 Fsub="AorN",
                                 Gnorm="AorN",
                                 Gsub="AorN",
				 Fnorm.var="AorN",
				 Fsub.var="AorN",
				 Gnorm.var="AorN",
				 Gsub.var="AorN"))

