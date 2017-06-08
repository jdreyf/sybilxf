#'@export convert_units
#'@title Convert pmol/min/ug_protein or pmol/min/ug_DNA to mmol/gDw/hr, the units of the metabolic model
#'@author Alfred Ramirez
#'@description This function takes the matrix returned by \code{\link{sample_seahorse}} and multiplies it by a factor
#'@param x A matrix
#'@param mwr Macromolecular weight ratio
#'@details The macromolecular weight ratio is the relationship between the initial normalization 
#'(typically DNA or protein) and total dry weight. It can vary substantially depending on cell type. By default,
#'mwr = 0.706, which is the estimated protein/dry_weight ratio from the biomass reaction in Recon 2. For adipocytes,
#'mwr = 0.026

convert_units <- function(x, mwr = 0.706){
  #Units are initially pmol/min/ug_macromolecule
  
  #Multiply by 1e-9 to convert from pmol to mmol
  c1 <- 1e-9
  
  #Multiply by 60 to convert from min to hrs
  c2 <- 60
  
  #Multiply by the mcr to convert from ug_macromolecule to ugDw
  c3 <- mwr
  
  #Multiply by 1e6 to from ugDw to gDw
  c4 <- 1e6
  
  y <- x*c1*c2*c3*c4
  y
}

