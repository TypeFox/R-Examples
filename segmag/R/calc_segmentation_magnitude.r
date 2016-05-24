calc_segmentation_magnitude <- function(segmag)
{
  # Baut ein Array in das die Segmentierungsstärke geschrieben wird
  # Um jeden Tastendruck wird Gauss gelegt
  # Gauss vorberechnet mit cutoff
  # Eine VPn jeweils Hüllfunktion, damit maximaler Beitrag beschränkt ist
  # Statt in Zeit wird in Indizes des Arrays gerechnet, da schneller geht.. Dazu alle mittels time_steps umgerechnet, also 1 / time_steps Arrayfelder je Sekunde

  if (! is.segmag(segmag)) stop("segmag must be an object of class segmag")
  
  # Vektor mit Segmentierungsstärke über Zeit in time_steps als Indizes
  segmentation_magnitude_overall <- numeric(segmag$index_time_max+1)
  
  for (id in levels(segmag$ids))
  {    
    index_keypresses <- segmag$index_keypresses[segmag$ids == id]
    
    calc_segmentation_magnitude_impl(segmentation_magnitude_overall,index_keypresses,segmag$gauss_values,segmag$gauss_n_indexes_per_side,segmag$indexes_gauss_offset)
  }
  
  return( data.frame(time=seq(segmag$time_min, segmag$time_min + (segmag$index_time_max*segmag$time_steps), segmag$time_steps), segmentation_magnitude=segmentation_magnitude_overall) )
}