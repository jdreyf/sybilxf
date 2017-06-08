#'@export fluxvar_rewrite
#'@title Rewrite fluxVar result into a matrix
#'@author Alfred Ramirez
#'@description This function rewrites an object of \code{\link[sybil]{optsol}} into a matrix
#'@param x A object of class \code{\link[sybil]{optsol}} produced by fluxVar

fluxvar_rewrite <- function(x){
  results_vector <- x@lp_obj
  reactions_tested <- x@react@react_id
  
  output_mat <- matrix(results_vector, byrow=F, ncol=2, nrow=length(reactions_tested))
  colnames(output_mat) <- c("Min", "Max")
  rownames(output_mat) <- reactions_tested
  output_mat
}
