#'@title Map the seahorse fluxes to metabolic reactions
#'@author Alfred Ramirez, Jonathan Dreyfuss
#'@description This function takes the matrix of sampled seahorse measurements returned by \code{\link{sample_seahorse}}
#'and maps the fluxes to metabolic reactions. It returns a matrix with the mapped fluxes of
#'basal oxygen consumption, mitochondrial oxygen consumption, mitochondrial ATP production,
#'mitochondrial ATP leak, and basal extracellular acidification
#'@param x A matrix of sampled seahorse measurements
#'@param model.nm The metabolic model name. One of "2.1A", "2.1x", or "2.2".
#'@export

map_seahorse <- function(x, model.nm=c("2.1A", "2.1x", "2.2")){
  stopifnot(model.nm %in% c("2.1A", "2.1x", "2.2"))
  # These reaction names are specific to Recon 2.1A & Recon 2.1x.
  if (model.nm %in% c("2.1A", "2.1x")){
    exp_coefs <- c(o2_in="EX_o2(e)in", o2_out="EX_o2(e)ex", atp_synth="ATPS4m", atp_demand="DM_atp_m_", 
                   o2_trans_mit="O2tm", h_in="EX_h(e)in", h_out="EX_h(e)ex")
  }
  if (model.nm=="2.2"){
    exp_coefs <- c(o2_ex="EX_o2(e)", atp_synth="ATPS4m", atp_demand="DM_atp_c_", 
                   o2_trans_mit="O2tm", h_ex="EX_h(e)")
  }
  
  ub <- setNames(paste0(exp_coefs,"_ub"), nm=names(exp_coefs))
  lb <- setNames(paste0(exp_coefs,"_lb"), nm=names(exp_coefs))
  
  output_mat <- matrix(0, ncol=ncol(x), nrow=length(exp_coefs)*2)
  colnames(output_mat) <- colnames(x)
  rownames(output_mat) <- c(ub, lb)
  
  #set some fluxes to a value by setting both their upper and lower bounds to that value
  if (model.nm=="2.2"){
    output_mat[lb["o2_ex"],] <- -1*x["OCR_basal",]
    output_mat[ub["o2_ex"],] <- -1*x["OCR_basal",]
  } else {
    output_mat[lb["o2_in"],] <- x["OCR_basal",]
    output_mat[ub["o2_in"],] <- x["OCR_basal",]
    output_mat[lb["o2_out"],] <- 0
    output_mat[ub["o2_out"],] <- 0
  }
  output_mat[lb["atp_synth"],] <- 4.6*(x["OCR_basal",] - x["OCR_oligo",])
  output_mat[ub["atp_synth"],] <- 4.6*(x["OCR_basal",] - x["OCR_oligo",])
  output_mat[lb["atp_demand"],] <- 4.6*(x["OCR_oligo",] - x["OCR_rotenone",])
  output_mat[ub["atp_demand"],] <- 4.6*(x["OCR_oligo",] - x["OCR_rotenone",])
  output_mat[lb["o2_trans_mit"],] <- x["OCR_basal",] - x["OCR_rotenone",]
  output_mat[ub["o2_trans_mit"],] <- x["OCR_fccp",] - x["OCR_rotenone",]
  if (any(output_mat[ub["o2_trans_mit"],] < output_mat[lb["o2_trans_mit"],])){
    warning("Some OCR_fccp sampled fluxes were less than OCR_basal sampled fluxes. To avoid infeasibility, 
            these were set to the OCR_basal values.")
    wh.cols <- which(output_mat[ub["o2_trans_mit"],] < output_mat[lb["o2_trans_mit"],])
    output_mat[ub["o2_trans_mit"], wh.cols] <- output_mat[lb["o2_trans_mit"], wh.cols]
  }
  
  #Lactate and extracellular acid may be negative and the associated reactions have different 
  #stoichiometries due to the total carbon constraint, thus an if statement is needed to check to determine which one to set.
  #Ideally, there should never be a case where the PPR is negative (except a few edges cases in biology).
  if (model.nm=="2.2"){
    for(i in 1:ncol(x)){
      output_mat[lb["h_ex"],i] <- x["PPR_basal",i]
      output_mat[ub["h_ex"],i] <- x["PPR_basal",i]
    }
  } else {
    for(i in 1:ncol(x)){
      if(x["PPR_basal",i] < 0){
        output_mat[lb["h_in"],i] <- abs(x["PPR_basal",i])
        output_mat[ub["h_in"],i] <- abs(x["PPR_basal",i])
        output_mat[lb["h_out"],i] <- 0
        output_mat[ub["h_out"],i] <- 0
      } else {
        output_mat[lb["h_in"],i] <- 0
        output_mat[ub["h_in"],i] <- 0
        output_mat[lb["h_out"],i] <- x["PPR_basal",i]
        output_mat[ub["h_out"],i] <- x["PPR_basal",i]
      }
    }
  }#end else model not reversible

  #if use irrev model, then no fluxes should be negative
  if(model.nm %in% c("2.1A", "2.1x") && any(output_mat < 0 )){
    warning("Model is irreversible so fluxes are expected to be positive, but some sampled fluxes were negative.")
  }
  output_mat
}