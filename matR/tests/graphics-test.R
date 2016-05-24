library(matR)

#-----------------------------------------------------------------------------------------
#  OK FOR CRAN
#-----------------------------------------------------------------------------------------

parMap (xx1, 											# extraneous element in map
		name.map = c (									# no map for "col"
			col="host_common_name",						# full map for "pch" includes "NA"
			pch="samp_store_temp"),
		value.map = list (
			pch = c (
				"-80"="+",
				"NA"="x"),
			extraneous = 5000))
parMap (xx1, 											# extraneous element in map
		name.map = c (									# partial map for "col" (1 specified)
			col="host_common_name",						# full map for "pch" includes "NA"
			pch="samp_store_temp"),
		value.map = list (
			col = c (Mouse="blue"),
			pch = c (
				"-80"="+",
				"NA"="x"),
			extraneous = 5000))
parMap (xx1, 											# extraneous element in map
		name.map = c (									# partial map for "col" (2 specified)
			col="host_common_name",						# full map for "pch" includes "NA"
			pch="samp_store_temp"),
		value.map = list (
			col = c (
				Mouse="blue",
				cow="red"),
			pch = c (
				"-80"="+",
				"NA"="x"),
			extraneous = 5000))
parMap (xx1, 											# extraneous element in map
		name.map = c (									# partial map for "col", "NA" spec for "other"
			col="host_common_name",						# full map for "pch" includes "NA"
			pch="samp_store_temp"),
		value.map = list (
			col = c (
				Mouse="blue",
				"NA"="brown"),
			pch = c (
				"-80"="+",
				"NA"="x"),
			extraneous = 5000))
