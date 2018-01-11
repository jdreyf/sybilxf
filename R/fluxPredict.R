#'@title Make predictions from seahorse data by sampling and minimizing total flux
#'@author Alfred Ramirez, Jonathan Dreyfuss
#'@description This function integrates the sampled seahorse measurements as constraints into the specified model, 
#'optionally maximizes an objective reaction, minimizes total flux for each sample, and returns a matrix 
#'of reactions-by-samples where the entries are predicted fluxes.
#'@param model An object of class \code{\link[sybil]{modelorg}}
#'@param seahorse_data A data.frame returned by \code{\link{map_seahorse}}
#'@param biomass_est Estimated biomass flux. Default = 0.
#'@param alg Either "fba" (default) to optimize a reaction or "mtf" to only minimize total flux (and not do fba).
#'@param obj_rxn Reaction in \code{model} to optimize if \code{alg="fba"}. Default \code{NULL} optimizes ATP demand reaction.
#'@param model.nm The metabolic model name. One of "2.1A", "2.1x", or "2.2".
#'@details For parallel computing, a parallel backend must be registered. See \code{\link[foreach]{foreach}} for details.
#'@export

fluxPredict <-function(model, seahorse_data, biomass_est=0, alg=c("fba", "mtf"), 
                       obj_rxn=NULL, model.nm, solver=c("gurobi", "glpk")){
  alg <- match.arg(alg)
  solver <- match.arg(solver)
  stopifnot(model.nm %in% c("2.1A", "2.1x", "2.2"))
  if (alg=="fba" && !is.null(obj_rxn) && !(obj_rxn %in% sybil::react_id(model))){
    stop("obj_rxn ", obj_rxn, " not in react_id(model).")
  }
  if (alg=="mtf" && !is.null(obj_rxn)) warning("alg=mtf, so obj_rxn will be ignored.")
  if (!is.numeric(biomass_est)||biomass_est<0) stop("biomass must be numeric and non-negative.")
  
  #biomass reaction has same name in 2.1A & 2.2
  biomass <- "biomass_reaction"
  #exp_coefs <- c("EX_o2(e)in","EX_o2(e)ex","ATPS4m", "DM_atp_m_", "O2tm", "EX_lac_L(e)in", "EX_lac_L(e)ex", "biomass_reaction")
  exp_coefs <- c(gsub("_lb$", "", grep("_lb$", rownames(seahorse_data), value=TRUE)), biomass=biomass)
  #ub and lb should match
  exp_coefs_ub <- c(gsub("_ub$", "", grep("_ub$", rownames(seahorse_data), value=TRUE)), biomass=biomass)
  if (any(exp_coefs!=exp_coefs_ub)) stop("seahorse_data should have matching upper bounds and lower bounds, but it does not.")
  
  #set obj
  if (alg=="fba" && is.null(obj_rxn)){
    #atp demand reaction has this form in 2.1A & 2.2
    dm_atp <- grep("DM_atp", exp_coefs, value=TRUE)
    if (length(dm_atp)>1){
      stop("Expected one ATP demand reaction with 'DM_atp' in name from rownames(seahorse_data), but found several.")
    }
    obj_rxn <- dm_atp
  }#end if fba & obj_rxn
  
  #create wts to nudge to co2 + lactate efflux
  if (alg=="mtf"){
    wts <- rep(1, times=length(model@react_id))
    names(wts) <- model@react_id
    #found rxns initially in vmh
    #then found rxns in 2.1A by using grep w/ fixed=TRUE
    #L_LACt2r transports lac c <-> e
    if (model.nm=="2.2"){
      lac.rxns <- c("L_LACt2r", "EX_lac_L(e)", "LDH_L")
      co2.rxns <- c("CO2t", "CO2tm", "EX_hco3(e)", "H2CO3D", "H2CO3Dm", "r0941", "r1418")
      model <- sybil::changeBounds(model, react=c("H2CO3D", "H2CO3Dm", "EX_lac_L(e)", "EX_hco3(e)"), lb=0)
      model <- sybil::changeBounds(model, react=c("L_LACt2r", "LDH_L"), ub=0)
      #no co2 ex
      model <- sybil::changeBounds(model, react="EX_co2(e)", lb=0, ub=0)
      
    } else if (model.nm %in% c("2.1A", "2.1x")){
      lac.rxns <- c("L_LACt2r", "EX_lac_L(e)ex", "LDH_L")
      co2.rxns <- c("CO2t", "CO2tm", "EX_hco3(e)ex", "H2CO3D", "r0941", "r1418")
    }
    low.wt.rxns <- c(lac.rxns, co2.rxns)
    stopifnot(low.wt.rxns %in% model@react_id)
    wts[low.wt.rxns] <- 0.01
  }#end create wts
  
  final_output <- foreach(i=1:ncol(seahorse_data), .combine=cbind) %dopar% {
    model_lb <- c(seahorse_data[grep("_lb$", rownames(seahorse_data)),i], biomass_est)
    model_ub <- c(seahorse_data[grep("_ub$", rownames(seahorse_data)),i], biomass_est)
    stopifnot(exp_coefs %in% model@react_id)
    model <- sybil::changeBounds(model, react=exp_coefs, lb=model_lb, ub=model_ub)
    
    if (alg=="fba"){
      model <- sybil::changeObjFunc(model, obj_rxn, obj_coef=1)
      model_test_flux <- sybil::optimizeProb(model, algorithm="fba", lpdir="max")  
      lpsolution <- model_test_flux@lp_stat 
      #If the constraints are not feasible, we return a column of NAs for those constraints
      if (lpsolution == 0){
        output <- NA
      } else {
        output_fluxes <- sybil::optimizeProb(model, alg="mtf", mtfobj=mod_obj(model_test_flux))
        output <- sybil::getFluxDist(output_fluxes)
      }
      output
    } else { 
      #if not fba, then min weighted flux to nudge to hco3 + lactate efflux
      # output_fluxes <- sybil::optimizeProb(model, alg="mtf")
      #output_fluxes <- sybil::optimizeProb(model, alg="fba", react=1:length(wts), obj_coef=wts)
      #output <- sybil::getFluxDist(output_fluxes)
      nrxns <- model@react_num
      nmets <- model@met_num
      eye <- simple_triplet_diag_matrix(v=rep(1, times=nrxns))
      neg.eye <- simple_triplet_diag_matrix(v=rep(-1, times=nrxns))
      mat1 <- abind_simple_sparse_array(as.simple_triplet_matrix(model@S), simple_triplet_zero_matrix(nrow = nmets, ncol=nrxns), MARGIN = 2)
      mat2 <- abind_simple_sparse_array(eye, neg.eye, MARGIN = 2)
      mat3 <- abind_simple_sparse_array(neg.eye, neg.eye, MARGIN = 2)                           
      mat <- abind_simple_sparse_array(mat1, mat2, mat3, MARGIN=1)
      mat <- as.simple_triplet_matrix(mat)
      
      obj <- c(numeric(nrxns), wts)
      rhs <- numeric(nmets+2*nrxns)
      dirs <- rep(c("==", "<="), times=c(nmets, 2*nrxns))
      if (solver=="glpk"){
        lb <- list(ind = 1:nrxns, val = model@lowbnd)
        ub <- list(ind = 1:nrxns, val = model@uppbnd)
        opt <- Rglpk::Rglpk_solve_LP(obj=obj, mat=mat, dir=dirs, rhs=rhs, bounds=list(lower=lb, upper=ub))
        if (opt$status==0) print("OPTIMAL SOLUTION FOUND")
        output <- opt$solution[1:nrxns]
      } else if (solver=="gurobi"){
        sense <- sub("==", "=", dirs)
        lb <- c(model@lowbnd, numeric(nrxns))
        ub <- c(model@uppbnd, rep(1000, nrxns))
        mod <- list(obj=obj, A=mat, rhs=rhs, sense=sense, lb=lb, ub=ub)
        opt <- gurobi::gurobi(model=mod)
        if (opt$status=="OPTIMAL") print("OPTIMAL SOLUTION FOUND")
        output <- opt$x[1:nrxns]
      }
    }
  }#end foreach
  rownames(final_output) <- sybil::react_id(model)
  final_output
}