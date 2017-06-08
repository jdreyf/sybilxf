#'@title Get all the reactions that contain a specified metabolite
#'@author Alfred Ramirez
#'@description This function returns the stoichiometric matrix of all the reactions containing the specified metabolite(s)
#'@param metabolites A character vector of metabolites
#'@param model A object of class \code{\link[sybil]{modelorg}}

get_reactions <- function(metabolites, model){
  s_mat <- S(model)
  colnames(s_mat) <- react_id(model)
  rownames(s_mat) <- met_id(model)
  
  search_space <- glob2rx(paste0(metabolites,"[*]"))
  mets <- rownames(s_mat)[grep(search_space, rownames(s_mat))]
  message("Found these metabolites corresponding to the input:", mets)
  s_mat <- s_mat[rownames(s_mat) %in% mets,]
  s_mat <- s_mat[,colSums(abs(s_mat)) > 0 ]
  s_mat
}
