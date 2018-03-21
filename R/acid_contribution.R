#' @title Calculate contributions of respiration and glycolysis to extracellular acid production
#' @author Jonathan Dreyfuss
#' @description Calculate contributions of respiration and glycolysis to extracellular acid production from 
#' Mookerjee et al. (2015a), equation 3.
#' @param PPR_total Proton Production Rate = ECAR / buffering power, from Mookerjee et al. (2015b), equation 3.
#' @param OCR_total Oxygen Consumption Rate (pmol O2/min).
#' @param OCR_nonmit Non-mitochondrial OCR remaining after complete mitochondrial inhibition, eg with rotenone.
#' @param mho maximum H+/O2. 1.0 for glucose and glycogen, 0.8 for pyruvate, 0.65 for palmitate.
#' @param pH pH of the experiment.
#' @param pK1 Overall pK for CO2(aq) + H2O -> HCO3− + H+, which is 6.093 at 37 Celsius.
#'
#' @references Mookerjee SA, Goncalves RL, Gerencser AA, Nicholls DG, Brand MD. The contributions of respiration and 
#' glycolysis to extracellular acid production. Biochim Biophys Acta. 2015 Feb;1847(2):171-81.
#' @references Mookerjee SA, Brand MD. Measurement and Analysis of Extracellular Acid Production to Determine 
#' Glycolytic Rate. J Vis Exp. 2015 Dec 12;(106):e53464.

acid_contribution <- function(PPR_total, OCR_total, OCR_nonmit, mho=1, pH=7.4, pK1=6.093){
  PPR_resp <- (OCR_total - OCR_nonmit)*mho*10^(pH-pK1)/(1 + 10^(pH-pK1))
  PPR_glyc <- PPR_total - PPR_resp
  return(c(PPR_resp=PPR_resp, PPR_glyc=PPR_glyc))
}