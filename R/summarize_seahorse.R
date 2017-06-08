#'@export summarize_seahorse
#'@title Summarize seahorse data by sample and region
#'@author Alfred Ramirez
#'@description This function takes the long format table exported by the XF Wave software for the 
#'oxidiative stress test and returns a matrix of means and standard deviations for each region
#'of the assay (namely basal, oligomycin, fccp, and rotenone).  The GroupName column must have the sample names.
#'The final measurement of each region is used for the mean and standard deviation across wells.
#'@param x A data.frame
#'@param n The number of regions

summarize_seahorse <- function(x, n=4){
  expected_colnames <- c("Measurement", "GroupName", "Time", "OCR", "OCR.Error", "ECAR", "ECAR.Error", "PPR", "PPR.Error")
  if(!all(colnames(x) %in% expected_colnames)){
    stop("Unrecognized column names.  Colnames must be as follows:",expected_colnames)
  }
  samples <- unique(x[,"GroupName"])
  
  # I use length(unique(foo)) rather than max(foo) since there may be a time when
  # a whole region should be excluded from the analysis creating a discrepancy b/w the two
  m <- length(unique((x[,"Measurement"])))/n #The number measurements in each of the n regions
  output_mat <- matrix(0,nrow=length(samples), ncol=(4*n))
  rownames(output_mat) <- samples
  
  #Kudos to those that can wrap their heads around this.
  for(i in 1:length(samples)){
    mat <- x[x[,"GroupName"] == samples[i],]
    for(j in 1:n){      
      mini_mat <- tail(mat[((j-1)*m+1):(j*m),],1)
      output_mat[i,(1:4)+4*(j-1)] <- as.matrix(mini_mat[,c("OCR", "OCR.Error", "PPR", "PPR.Error")])
    }
  }
    
  if(ncol(output_mat) == 16){  
    #Most people probably use oligo, fccp, rotenone
    colnames(output_mat) <- c("OCR_basal","OCR_sd_basal", "PPR_basal", "PPR_sd_basal",
                              "OCR_oligo", "OCR_sd_oligo", "PPR_oligo", "PPR_sd_oligo",
                              "OCR_fccp","OCR_sd_fccp", "PPR_fccp", "PPR_sd_fccp", 
                              "OCR_rotenone", "OCR_sd_rotenone", "PPR_rotenone", "PPR_sd_rotenone")
  } else if(ncol(output_mat) == 20){
    #We used pyruvate, oligo, fccp, rotenone
    colnames(output_mat) <- c("OCR_basal","OCR_sd_basal", "PPR_basal", "PPR_sd_basal",
                              "OCR_pyr", "OCR_sd_pyr","PPR_pyr","PPR_sd",
                              "OCR_oligo", "OCR_sd_oligo", "PPR_oligo", "PPR_sd_oligo",
                              "OCR_fccp","OCR_sd_fccp", "PPR_fccp", "PPR_sd_fccp", 
                              "OCR_rotenone", "OCR_sd_rotenone", "PPR_rotenone", "PPR_sd_rotenone")
  } else {
    #For everyone else
    colnames(output_mat) <- paste0(rep(c("OCR_region", "OCR_sd_region", "PPR_region", "PPR_sd_region"),times=n),rep(1:n,each=4))
  }
  output_mat
}
