
################################################################################
################## FICHIER TEST DE plspm.formula.R #############################
################################################################################

library(plspm.formula)


###############################################################################
################# TESTS WITH SATISFACTION DATASET #############################
###############################################################################

# Load dataset
data("plspmsat")
# Model specification by formulas 
satmodele <- "
            # measure model specification
              EXPE =~ expe1+expe2+expe3+expe4+expe5
              IMAG =~ imag1+imag2+imag3+imag4+imag5
              LOY =~ loy1+loy2+loy3+loy4
              SAT =~ sat1+sat2+sat3+sat4
              VAL =~ val1+val2+val3+val4 
              QUAL =~ qual1+qual2+qual3+qual4+qual5 
            # outer model specification 
              EXPE ~~ IMAG
              LOY ~~ IMAG+SAT
              SAT ~~ IMAG+EXPE+QUAL+VAL
              VAL ~~ EXPE+QUAL
              QUAL ~~ EXPE
          "
# PLSPM estimation based on formula
satmodes <- rep("A",6)  
satres.plspm <- plspm.formula(Formula = satmodele, Data = plspmsat, modes = satmodes, 
                             plot.outer = TRUE, plot.inner = TRUE, scaled = FALSE)
# computer the PLSPM parameters
sat.param <- plspm.params(Formula = satmodele, Data = plspmsat)
sat.param$inner.mat  # inner matrix
sat.param$outer.list  #  outer list



###############################################################################
################ TESTS WITH ADULTS RESILIENCE DATASET #########################
###############################################################################

# Load dataset
data(bkadulte)
# Model specification by formulas 
bkmodele <- "
             # modele de mesure 
              SDH =~ senshum+creativite 
              APS =~ communic+sociabilite+altruiste+relation
              HRP =~ planific+solution+autonome
              SCI =~ estime+confiance+favenir
              SPI =~ optimisme+persever+religion 
            # interactions 
              HRP ~~ SDH+SCI
              SCI ~~ SDH
              APS ~~ SDH
              SPI ~~ APS+SCI
          "
# PLSPM estimation based on formula
bkmodes <- rep("A",5)  
bkres.plspm <- plspm.formula(Formula = bkmodele, Data = bkadulte, modes = bkmodes, 
                             plot.outer = TRUE, plot.inner = TRUE, scaled = FALSE)
# Computation plspm parameters only based on formula
bkres.param <- plspm.params(Formula = bkmodele, Data = bkadulte)
bkres.param$inner.mat
bkres.param$outer.list
