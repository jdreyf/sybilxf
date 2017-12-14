#'@title Sample from the seahorse data
#'@author Alfred Ramirez
#'@description This function takes the data.frame returned by the function \code{\link{summarize_seahorse}} and returns
#'a data.frame with sampled measurements for the specified sample
#'@param x A data.frame
#'@param sample.nm The sample name 
#'@param nsamples The number of samples to generate
#'@export

sample_seahorse <- function(x, sample.nm, nsamples=150){
  if (!(sample.nm %in% rownames(x))){
    stop(paste0("The sample.nm '", sample.nm, "' must be in the rownames of x"))
  }
  
  output_mat <- matrix(0, ncol=nsamples, nrow=8)
  colnames(output_mat) <- paste0("sample", 1:nsamples)
  rownames(output_mat) <- c("OCR_basal", "OCR_oligo", "OCR_fccp", "OCR_rotenone", "PPR_basal", "PPR_oligo", "PPR_fccp", "PPR_rotenone")
  
  i <- sample.nm
  output_mat["OCR_basal",] <- rnorm(nsamples,mean=x[i,"OCR_basal"],sd=x[i,"OCR_sd_basal"])
  output_mat["OCR_oligo",] <- rnorm(nsamples,mean=x[i,"OCR_oligo"],sd=x[i,"OCR_sd_oligo"])
  output_mat["OCR_fccp",] <- rnorm(nsamples,mean=x[i,"OCR_fccp"],sd=x[i,"OCR_sd_fccp"])
  output_mat["OCR_rotenone",] <- rnorm(nsamples,mean=x[i,"OCR_rotenone"],sd=x[i,"OCR_sd_rotenone"])
  
  output_mat["PPR_basal",] <- rnorm(nsamples,mean=x[i,"PPR_basal"],sd=x[i,"PPR_sd_basal"])
  output_mat["PPR_oligo",] <- rnorm(nsamples,mean=x[i,"PPR_oligo"],sd=x[i,"PPR_sd_oligo"])
  output_mat["PPR_fccp",] <- rnorm(nsamples,mean=x[i,"PPR_fccp"],sd=x[i,"PPR_sd_fccp"])
  output_mat["PPR_rotenone",] <- rnorm(nsamples,mean=x[i,"PPR_rotenone"],sd=x[i,"PPR_sd_rotenone"])
  
  output_mat
}
