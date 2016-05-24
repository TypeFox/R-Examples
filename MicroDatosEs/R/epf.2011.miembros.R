###################################################################
# Diego Paniagua Sánchez
# 20150906
# Reads the members EPF microdata in its 2011 version
###################################################################

epf.2011.miembros <- function(file) 
{       
        file.column <- create.spss.column(system.file("metadata", 
                                                      "epf_2011_miembros_mdat1.txt", package = "MicroDatosEs"), system.file("metadata", 
                                                                                                                     "epf_2011_miembros_mdat2.txt", package = "MicroDatosEs"), fileEncoding = "UTF-8")
        
        
        file.var <- create.spss.var(system.file("metadata", "epf_2011_miembros_mdat1.txt", 
                                                package = "MicroDatosEs"), fileEncoding = "UTF-8")
        file.vals <- create.spss.vals(system.file("metadata", "epf_2011_miembros_mdat2.txt", 
                                                  package = "MicroDatosEs"), fileEncoding = "UTF-8")
        file.missing <- system.file("metadata", "epf_2011_miembros_mdat3.txt", 
                                    package = "MicroDatosEs")
        x <- spss.fixed.file(file = file, columns.file = file.column, 
                             varlab.file = file.var, missval.file = file.missing, 
                             codes.file = file.vals)
        
        as.data.set(x)
}

