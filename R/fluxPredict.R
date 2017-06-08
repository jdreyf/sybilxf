#'@export fluxPredict
#'@import "foreach" "sybil"
#'@title Make predictions from seahorse data by sampling and minimizing total flux
#'@author Alfred Ramirez
#'@description This function integrates the sampled seahorse measurements as constraints into the specified model,
#'minimizes total flux for each sample, and returns a matrix of reactions x samples where the entries are predicted fluxes.
#'@param model An object of class \code{\link[sybil]{modelorg}}
#'@param seahorse_data A data.frame returned by \code{\link{map_seahorse}}
#'@param biomass_est Estimated biomass flux.  Not required, by may be useful.
#'@param cores The number of cores to use for parallel computing.  If cores is greater than 1, a parallel backend
#'must be registered.

fluxPredict <-function(model, seahorse_data, biomass_est=0, cores=1){
  final_output <- foreach(i=1:ncol(seahorse_data), .combine=cbind) %dopar% {
    exp_coefs <- c("EX_o2(e)in","EX_o2(e)ex","ATPS4m", "DM_atp_m_", "O2tm", "EX_lac_L(e)in", "EX_lac_L(e)ex", "biomass_reaction")  
    model_lb <- c(seahorse_data[grep("_lb$", rownames(seahorse_data)),i], biomass_est)
    model_ub <- c(seahorse_data[grep("_ub$", rownames(seahorse_data)),i], biomass_est)
    model <- changeBounds(model, react= exp_coefs, lb=model_lb, ub=model_ub)
    model <- changeObjFunc(model, "DM_atp_m_", obj_coef=1)    
    
    #If the constraints are not feasible, we return a column of NAs for those constraints
    #for those constraints
    
    model_test_flux <- optimizeProb(model, algorithm="fba", lpdir="max")  
    lpsolution <- model_test_flux@lp_stat 
    
    if(lpsolution == 0){
      output <- NA
    } else{
      output_fluxes <- optimizeProb(model, alg="mtf", mtfobj=mod_obj(model_test_flux))
      output <- getFluxDist(output_fluxes)
    }
    output
  }
  rownames(final_output) <- react_id(model)
  final_output
}