FPKMToTPM <- function() {
  
  # Process one column at a time.
  tpm <-  {
    counts <- data.frame(fread("GSE118612_counts.txt", sep = "\t"), check.names = F, row.names = 1)
    rate = sum(counts)
    (counts/rate) * 1e6
  }
  
  log2TPM = log2(tpm+1)
  # Copy the row and column names from the original matrix.
  colnames(log2TPM) <- colnames(counts)
  rownames(log2TPM) <- rownames(counts)
  
  fwrite(log2TPM, "GSE61220_EMTtpm.txt", sep = "\t", row.names = T)
}

FPKMToTPM()

