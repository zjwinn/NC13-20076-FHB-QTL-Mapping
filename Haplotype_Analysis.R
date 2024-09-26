# Read in study data
load("Final_Analysis.RData")

# Lines 
lines <- as.character(cross_file$pheno$GENOTYPE)

# Pull marker matrix
vcf <- gaston::read.vcf("NC13-20076xGA06493-13LE6_filt_imp.vcf.gz", convert.chr = FALSE)

# Add parents 
lines <- c(lines, "13955-GA06493-13LE6", "13955-NC13-20076", "13955-NCAG11", "13955-AGS-2026", "13955-JAMESTOWN") 

# Pull lines
vcf <- vcf[vcf@ped$id %in% lines,]

# Make matrix
vcf <- gaston::as.matrix(vcf)

# Make list
qtl <- "QDon.nc-1A"

# For i in QTL
haplotype <- function(i, qtl_table, vcf, threshold = 0.99){
  # Pull Data
  chr <- qtl_table[qtl_table$QTL==i,"Chromosome"]
  left_pos <- qtl_table[qtl_table$QTL==i,"Left Position (Mbp)"] * 1000000
  right_pos <- qtl_table[qtl_table$QTL==i,"Right Position (Mbp)"] * 1000000
  
  # Pull marker matrix
  marker_info <- data.frame(marker = colnames(vcf),
                            chromosome = substr(colnames(vcf),2,3),
                            position = as.numeric(gsub(".*_", "", colnames(vcf))))
  marker_info <- marker_info[marker_info$chromosome==chr & 
                               marker_info$position>left_pos &
                               marker_info$position<right_pos,]
  marker_matrix <- as.matrix(vcf[,colnames(vcf) %in% marker_info$marker])
  
  # Correlations
  corr <- cor(marker_matrix)
  corr <- caret::findCorrelation(abs(corr), cutoff = threshold)
  
  # Remove those high LD markers
  marker_matrix <- marker_matrix[,-corr]
  
  # do PCA
  pca <- summary(prcomp(t(marker_matrix)))
  select_pca <- as.data.frame(t(pca$importance))
  select_pca <- select_pca[select_pca$`Proportion of Variance`>=.01,]
  pca <- pca$x
  
  # do nbclust
  nbc <- suppressWarnings(NbClust::NbClust(data = pca[,c(rownames(select_pca))], distance = "euclidean", method = "ward.D2"))
  
  # Pull proper clustering of markers
  proper_clustering <- nbc$Best.partition
  
  # Cluster results
  classification <- data.frame(Line = rownames(marker_matrix), check.names = FALSE, row.names = NULL)
  
  # Haploblocks
  haploblocks <- unique(proper_clustering)
  haploblocks <- haploblocks[order(haploblocks)]
  
  for(j in haploblocks){
    # Select markers in proper clustering
    temp1 <- marker_matrix[,names(proper_clustering[proper_clustering==j])]
    
    # do PCA
    temp2 <- summary(prcomp(temp1))
    temp3 <- as.data.frame(t(temp2$importance))
    temp3 <- temp3[temp3$`Proportion of Variance`>=.01,]
    temp2 <- temp2$x
    
    # do nbclust
    temp4 <- suppressWarnings(NbClust::NbClust(data = temp2[,c(rownames(temp3))], distance = "euclidean", method = "ward.D2"))
    
    # Pull proper clustering of markers
    temp5 <- data.frame(Line = names(temp4$Best.partition),
                        Cluster = temp4$Best.partition, 
                        check.names = FALSE,
                        row.names = NULL)
    
    # Rename column
    colnames(temp5)[2] <- paste("HaploBlock", j, sep = "_")
    
    # Bind in
    classification <- dplyr::left_join(classification, temp5, by = "Line")
    
    # Clean up
    remove(temp1, temp2, temp3, temp4, temp5)
  }
  
  # Return
  return(classification)
}

# Get haplotypes
haplotypes <- haplotype(qtl, qtl_table = vr_fdk_don_table_me, vcf = vcf, threshold = 0.99)

# Pull resistant and susceptible haplotypes
resistant <- haplotypes[haplotypes$Line=="13955-NC13-20076",]
resistant <- paste(resistant[,2:4], collapse = "_")
susceptible <- haplotypes[haplotypes$Line=="13955-GA06493-13LE6",]
susceptible <- paste(susceptible[,2:4], collapse = "_")

# Pull lines which are resistant and susceptible
haplotypes$Haplotype <- paste(haplotypes$HaploBlock_1, haplotypes$HaploBlock_2, haplotypes$HaploBlock_3, sep = "_")
haplotypes$Call <- ifelse(haplotypes$Haplotype %in% resistant, "R", ifelse(haplotypes$Haplotype %in% susceptible, "S", "UNKNOWN"))

# Make selection
selection <- haplotypes[haplotypes$Call!="UNKNOWN",]

# Regress to validate call
validation <- selection[,c("Line", "Call")]
colnames(validation)[1] <- "GENOTYPE"
validation <- dplyr::left_join(pheno, validation, by = "GENOTYPE")
validation <- tidyr::drop_na(validation, DON, Call)

# Regress
summary(lm(DON~Call, data = validation))

# Plot quick
plot(as.factor(validation$Call), validation$DON)

# Write out haplotypes
write.csv(haplotypes,
          file.path(getwd(), "results_csv_files/QDon.nc-1A_haplotypes.csv"),
          row.names = FALSE)
