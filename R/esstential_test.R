#'@title Block the specified reactions one at a time and maximiize the objective function
#'@author Alfred Ramirez
#'@description This function iterates through each specified reaction, sets the bounds to lb=ub=0, and 
#'maxmizes the objective function.
#'@param model A object of class \code{\link[sybil]{modelorg}}
#'@param reactions A character vector of reaction ids
#'@param solver The solver to use.  Default SYBIL_SETTINGS("SOLVER")

essential_test <- function(model, reactions, solver=SYBIL_SETTINGS("SOLVER")){
  output <- vector(mode="numeric", length=length(reactions))
  names(output) <- reactions
  for(i in 1:length(reactions)){
    print(paste0("Iterating through reaction:", reactions[i]))
    revised_model <- changeBounds(model, react=reactions[i], lb=0, ub=0)
    predicted_fluxes <- optimizeProb(revised_model, algorithm="fba", solver=solver)
    output[i] <- mod_obj(predicted_fluxes)
  }
  output
}
