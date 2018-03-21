#'@title Sample from the seahorse data
#'@author Alfred Ramirez, Jonathan Dreyfuss
#'@description This function takes the data frame returned by the function \code{\link{summarize_seahorse}} and returns
#'a data frame with sampled measurements for the specified sample
#'@param x A data.frame
#'@param sample.nm The sample name 
#'@param nsamples The number of samples to generate
#'@param nrep The number of replicates, used to calculate standard error
#'@export

sample_seahorse <- function(x, sample.nm, nsamples=150, nrep=1){
  if (!(sample.nm %in% rownames(x))){
    stop(paste0("The sample.nm '", sample.nm, "' must be in the rownames of x"))
  }
  
  injections <- sub("^OCR_sd_", "", grep("^OCR_sd", colnames(x), value=TRUE))
  n.inj <- length(injections)
  message("Injections inferred from colnames(x) are: ", paste(injections, collapse=", "))
  
  output_mat <- matrix(0, ncol=nsamples, nrow=n.inj*2)
  colnames(output_mat) <- paste0("sample", 1:nsamples)
  #rownames(output_mat) <- c("OCR_basal", "OCR_oligo", "OCR_fccp", "OCR_rotenone", "PPR_basal", "PPR_oligo", "PPR_fccp", "PPR_rotenone")
  rownames(output_mat) <- paste(rep(c("OCR", "PPR"), each=n.inj), injections, sep="_")
  
  for (inj.nm in injections){
    output_mat[paste0("OCR_", inj.nm),] <- rnorm(nsamples, mean=x[sample.nm, paste0("OCR_", inj.nm)], 
                                                sd=x[sample.nm, paste0("OCR_sd_", inj.nm)]/sqrt(nrep))
    output_mat[paste0("PPR_", inj.nm),] <- rnorm(nsamples, mean=x[sample.nm, paste0("PPR_", inj.nm)], 
                                                sd=x[sample.nm, paste0("PPR_sd_", inj.nm)]/sqrt(nrep))
  }
  
  output_mat
}
