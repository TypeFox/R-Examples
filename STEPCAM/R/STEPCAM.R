fallout.dispersal <- function(new_community, new_richness) # function with dispersal event. Removal chance inversely depends of frequency species.
{
  new_community_rows <- sample(c(1:nrow(new_community)),new_richness,prob=new_community[,ncol(new_community)],replace=F);
  new_community <- new_community[new_community_rows,];
  return(new_community)
}

fallout.filtering <- function(new_community, filtering_fallout, new_richness,optimum,n_traits) # function with filtering event.
{
  new_community <- new_community[,-ncol(new_community)];
  new_community_traits <- new_community[c(2:(1+n_traits))];
  traits_with_optimum <- (rbind(new_community_traits,optimum));
  distances_traits <- dist(traits_with_optimum);
  distances_traits2 <- as.matrix(distances_traits);
  distance_from_optimum <- distances_traits2[nrow(distances_traits2),];
  distances_ordered <- order(distance_from_optimum);
  end_trait <- new_richness + 1;
  new_community <- new_community[c(distances_ordered[2:end_trait]),];
  
  return(new_community)
}

fallout.competition <- function(new_community, competition_fallout,n_traits) # function with limiting similarity event
{
  for(i in 1:competition_fallout){
   new_community <- new_community[,-ncol(new_community)];
   trait_distances2 <- as.matrix(dist(new_community[,c(2:(n_traits+1))]))
   trait_distances <- dist(new_community[,(2:(n_traits+1))])
   which(trait_distances2==min(trait_distances), arr.ind=TRUE)->m3
   as.vector(m3[,1])->a
   mina1 <- min(trait_distances2[a[1],c(-a[2],-a[1])]) 
   mina2 <- min(trait_distances2[a[2],c(-a[2],-a[1])])   
   min_overall <- which(c(mina1,mina2)==min(mina1,mina2))  
   species_out <- a[min_overall]  
	 new_community <- new_community[c(-species_out),]
   }

 return(new_community)
}

## the Kraft.generator function: function that runs (hybrid) STEPCAMs
STEPCAM<- function(params, species, abundances,taxa,esppres,community_number,n_traits,species_fallout)
{

 output <- matrix(nrow=taxa) 
 allfinaloutput <- matrix(nrow=taxa) # matrix in which all output will be written
 new_traits <- c()

 dispersal_fallout <- params[1] # species falling out through dispersal events
 filtering_fallout <- params[2] # species falling out through filtering events
 competition_fallout <- params[3] # species falling out through limiting similarity events

 if(is.na(dispersal_fallout)) dispersal_fallout <- 0;  #just in case something went wrong
 if(is.na(competition_fallout)) competition_fallout <- 0;
 if(is.na(filtering_fallout)) filtering_fallout <- 0;

 ## traits and abundances in community ##
 tr <- species[esppres,];  

 ## calculated trait mean. This mean will be used as the 'optimal' trait value in a community for the filtering model. Species very dissimilar from this value will be filtered out ##
 optimum <- c()  
 for(i in 1:n_traits){
  optimum[i] <- mean(tr[,(i+1)])
 }

 speciesnames <- row.names(species)
 new_community <- species
 
 fallout <- c(rep(1,dispersal_fallout),rep(2,filtering_fallout),rep(3,competition_fallout)) # vector containing all community assembly steps: 1 = dispersal event, 2 = filtering event, 3 = limiting similarity event
 
 for(i in 1:species_fallout)
 {
   new_richness <- nrow(species) - i
   type = fallout[i];
   if(type == 1){
      new_community <- fallout.dispersal(new_community,new_richness);
   }
   if(type == 2){
      new_community <- fallout.filtering(new_community,1,new_richness,optimum,n_traits);
      freq <- species$freq[as.numeric(row.names(new_community))]
      new_community <- cbind(new_community,freq)
   }
   if(type == 3) {
      new_community <- fallout.competition(new_community,1,n_traits);
      freq <- species$freq[as.numeric(row.names(new_community))]
      new_community <- cbind(new_community,freq)
   }
 }
 species_names_left <- as.numeric(row.names(new_community))
 new_community <- cbind(new_community) ; names(new_community)[ncol(new_community)] <- "abund"
 species_presences <- rep(0,taxa) ; species_presences[new_community[,1]] <- 1
 species_presences <- matrix(species_presences,nrow=taxa,ncol=1)
 ## bind all output ##
 output <- cbind(output,species_presences)
 allfinaloutput <- cbind(allfinaloutput,output)
 allfinaloutput <- allfinaloutput[,c(-1,-2)]

 return(allfinaloutput) # final modelled community
}