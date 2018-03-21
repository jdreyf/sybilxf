#' @title Rewrite fluxVar result into a matrix
#' @author Alfred Ramirez
#' @description This function rewrites an object of \code{\link[sybil]{optsol}} into a matrix
#' @param x A object of class \code{\link[sybil]{optsol}} produced by fluxVar
#' @references Ramirez AK, Lynes MD, Shamsi F, Xue R, Tseng YH, Kahn CR, Kasif S, Dreyfuss JM. Integrating Extracellular 
#' Flux Measurements and Genome-Scale Modeling Reveals Differences between Brown and White Adipocytes. Cell Rep 2017 
#' Dec; 21(11): 3040-3048.
#' @export


fluxvar_rewrite <- function(x){
  results_vector <- x@lp_obj
  reactions_tested <- x@react@react_id
  
  output_mat <- matrix(results_vector, byrow=F, ncol=2, nrow=length(reactions_tested))
  colnames(output_mat) <- c("Min", "Max")
  rownames(output_mat) <- reactions_tested
  output_mat
}
