#'@export map_seahorse
#'@title Map the seahorse fluxes to metabolic reactions
#'@author Alfred Ramirez, Jonathan Dreyfuss
#'@description This function takes the matrix of sampled seahorse measurements returned by \code{\link{sample_seahorse}}
#'and maps the fluxes to metabolic reactions. It returns a matrix with the mapped fluxes of
#'basal oxygen consumption, mitochondrial oxygen consumption, mitochondrial ATP production,
#'mitochondrial ATP leak, and basal extracellular acidification
#'@param x A matrix

map_seahorse <- function(x){
  # These reaction names are specific to Recon 2.1x.  Probably won't work for Recon 2.2.
  exp_coefs <- c("EX_o2(e)in","EX_o2(e)ex","ATPS4m", "DM_atp_m_", "O2tm", "EX_lac_L(e)in", "EX_lac_L(e)ex")  
  exp_ub <- paste0(exp_coefs,"_ub")
  exp_lb <- paste0(exp_coefs,"_lb")
  
  output_mat <- matrix(0, ncol=ncol(x), nrow=length(exp_coefs)*2)
  colnames(output_mat) <- colnames(x)
  rownames(output_mat) <- c(exp_ub, exp_lb)
  
  #set some fluxes to a value by setting both their upper and lower bounds to that value
  output_mat["EX_o2(e)in_lb",] <- x["OCR_basal",]
  output_mat["EX_o2(e)in_ub",] <- x["OCR_basal",]
  output_mat["EX_o2(e)ex_lb",] <- 0
  output_mat["EX_o2(e)ex_lb",] <- 0
  output_mat["ATPS4m_lb",] <- 4.6*(x["OCR_basal",] - x["OCR_oligo",])
  output_mat["ATPS4m_ub",] <- 4.6*(x["OCR_basal",] - x["OCR_oligo",])
  output_mat["DM_atp_m__lb",] <- 4.6*(x["OCR_oligo",] - x["OCR_rotenone",])
  output_mat["DM_atp_m__ub",] <- 4.6*(x["OCR_oligo",] - x["OCR_rotenone",])
  output_mat["O2tm_lb",] <- x["OCR_basal",] - x["OCR_rotenone",]
  output_mat["O2tm_ub",] <- x["OCR_fccp",] - x["OCR_rotenone",]
  if (any(output_mat["O2tm_ub",] < output_mat["O2tm_lb",])){
    warning("Some OCR_fccp sampled fluxes were less than OCR_basal sampled fluxes. To avoid infeasibility, 
            these were set to the OCR_basal values.")
    wh.cols <- which(output_mat["O2tm_ub",] < output_mat["O2tm_lb",])
    output_mat["O2tm_ub", wh.cols] <- output_mat["O2tm_lb", wh.cols]
  }
  
  #Lactate and extracellular acid may be negative and the associated reactions have different 
  #stoichiometries due to the total carbon constraint, thus an if statement is needed to check to determine which one to set.
  #Ideally, there should never be a case where the PPR negative (except a few edges cases in biology).
  for(i in 1:ncol(x)){
    if(x["PPR_basal",i] < 0){
      output_mat["EX_lac_L(e)in_lb",i] <- abs(x["PPR_basal",i])
      output_mat["EX_lac_L(e)in_ub",i] <- abs(x["PPR_basal",i])
      output_mat["EX_lac_L(e)ex_lb",i] <- 0
      output_mat["EX_lac_L(e)ex_ub",i] <- 0
    } else {
      output_mat["EX_lac_L(e)in_lb",i] <- 0
      output_mat["EX_lac_L(e)in_ub",i] <- 0
      output_mat["EX_lac_L(e)ex_lb",i] <- x["PPR_basal",i]
      output_mat["EX_lac_L(e)ex_ub",i] <- x["PPR_basal",i]
    }
  }

  if(any(output_mat < 0 )){
    warning("Negative values were found. Some samples may have values inconsistent with biological expectation")
  }
  output_mat
}
