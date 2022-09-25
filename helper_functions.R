# -----------------------------------------------------------------------------------------------------------------------------
# LOADING PACKAGES
# -----------------------------------------------------------------------------------------------------------------------------

library(TCGAbiolinks)
library(plyr)
library(ggplot2)
library(limma)
library(EDASeq)


# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR EXTRACTING INFORMATION ON BATCH, PATIENT AND ID.
# -----------------------------------------------------------------------------------------------------------------------------

get_IDs <- function(data) {
  IDs <- strsplit(c(colnames(data)), "-")
  IDs <- ldply(IDs, rbind)
  colnames(IDs) <- c('project', 'tss','participant', 'sample', "portion", "plate", "center")
  cols <- c("project", "tss", "participant")
  IDs$patient <- apply(IDs[,cols],1,paste,collapse = "-" )
  barcode <- colnames(data)
  IDs <- cbind(IDs, barcode)
  condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample))
  condition  <- gsub("01+[[:alpha:]]", "cancer", condition)
  IDs$condition <- condition
  IDs$myorder  <- 1:nrow(IDs)
  return(IDs)
}

# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR IDENTIFICATION OF SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES - LIMMA
# -----------------------------------------------------------------------------------------------------------------------------


DE_limma <- function(my.contrast, my.data, my.design, coLFC, coFDR) {
  fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design), my.contrast))
  tt <- topTable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
  
  #up <- data.frame(rownames(tt[tt$logFC >= coLFC & tt$adj.P.Val < coFDR, ]))
  #down <- data.frame(rownames(tt[tt$logFC <= -coLFC & tt$adj.P.Val < coFDR, ]))
  
  #colnames(up) <- as.character("up")
  #colnames(down) <- as.character("down")
  
  #tt <- subset(tt,abs(tt$logFC) >= coLFC & tt$adj.P.Val < coFDR)
  index.up <- which(tt$logFC >= coLFC & tt$adj.P.Val < coFDR)
  index.down <- which(tt$logFC <= -coLFC & tt$adj.P.Val < coFDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)
  
  #final <- list(up, down)
  
  return(tt)
}


