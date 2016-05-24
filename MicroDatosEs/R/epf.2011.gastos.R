###################################################################
# Diego Paniagua Sánchez
# 20150906
# Reads the expenses EPF microdata in its 2011 version
###################################################################


epf.2011.gastos <- function(file) 
{
        file.column <- create.spss.column(system.file("metadata", "epf_2011_gastos_mdat1.txt", package = "MicroDatosEs"), 
						system.file("metadata", "epf_2011_gastos_mdat2.txt", package = "MicroDatosEs"), fileEncoding = "UTF-8")
        file.var <- create.spss.var(system.file("metadata", "epf_2011_gastos_mdat1.txt", package = "MicroDatosEs"), fileEncoding = "UTF-8")
        file.vals <- create.spss.vals(system.file("metadata", "epf_2011_gastos_mdat2.txt", package = "MicroDatosEs"), fileEncoding = "UTF-8")
        file.missing <- system.file("metadata", "epf_2011_gastos_mdat3.txt", package = "MicroDatosEs")
        gastos_2011 <- spss.fixed.file(file = file, columns.file = file.column, 
                                        varlab.file = file.var, missval.file = file.missing, 
                                        codes.file = file.vals)
        fix.char.items(as.data.set(gastos_2011))
        as.data.set(gastos_2011)
}
