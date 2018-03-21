#' @title Summarize seahorse data by sample and region
#' @author Alfred Ramirez, Jonathan Dreyfuss
#' @description This function takes the long format table exported by the XF Wave software for the 
#' Mitochondrial stress test and returns a matrix of means and standard deviations for each region
#' of the assay (e.g. basal, oligomycin, fccp, and rotenone).  The GroupName column must have the sample names.
#' The final measurement of each region is used for the mean and standard deviation across wells.
#' @param x A data.frame
#' @param injections The names of the injections
#' @references Ramirez AK, Lynes MD, Shamsi F, Xue R, Tseng YH, Kahn CR, Kasif S, Dreyfuss JM. Integrating Extracellular 
#' Flux Measurements and Genome-Scale Modeling Reveals Differences between Brown and White Adipocytes. Cell Rep 2017 
#' Dec; 21(11): 3040-3048.
#' @export


summarize_seahorse <- function(x, injections=c("basal", "oligomycin", "fccp", "rotenone")){
  n <- length(injections)
  expected_colnames <- c("Measurement", "GroupName", "Time", "OCR", "OCR.Error", "PPR", "PPR.Error")
  if (!all(expected_colnames %in% colnames(x))){
    # stop("Unrecognized column names.  Colnames must be as follows:",expected_colnames)
    stop("Missing necessary columns:", paste(setdiff(expected_colnames, colnames(x)), collapse=", "))
  }
  samples <- unique(x[,"GroupName"])
  
  # I use length(unique(foo)) rather than max(foo) since there may be a time when
  # a whole region should be excluded from the analysis creating a discrepancy b/w the two
  m <- length(unique((x[,"Measurement"])))/n #The number measurements in each of the n regions
  output_mat <- matrix(0, nrow=length(samples), ncol=(4*n))
  rownames(output_mat) <- samples
  
  #Kudos to those that can wrap their heads around this.
  for(i in 1:length(samples)){
    mat <- x[x[,"GroupName"] == samples[i],]
    for (j in 1:n){      
      mini_mat <- tail(mat[((j-1)*m+1):(j*m),],1)
      output_mat[i,(1:4)+4*(j-1)] <- as.matrix(mini_mat[,c("OCR", "OCR.Error", "PPR", "PPR.Error")])
    }
  }
  
  colnames(output_mat) <- paste(c("OCR", "OCR_sd", "PPR", "PPR_sd"), rep(injections, each=4), sep="_")
  output_mat
}
