football <- heliumFootballs      # give it a shorter name
head(football,3)
 
football$diff <- football$Helium - football$Air        
t.test(football$diff)
