# Set the path to your directory
#data_dir <- "/Users/pathumnawarathna/Documents/Research/Simulation/keele/Probability_matrices_from_UNC/Genotype"
set_chr <- c(1:19,"X")

ccreate_additve_data <- function(data_dir, set_chr){
  options(digits = 15)
  full_dir <- paste0(data_dir,"/full")
  
  for (i in set_chr) {
    chromosome <- paste0("chr",i)
    
    # List all .RData files
    rdata_files <- list.files(path = paste0(full_dir,"/",chromosome,"/data/"), pattern = "\\.RData$", full.names = TRUE)
    
    # Load each file
    for (f in rdata_files) {
      
      file_name <- load(f)
      full_matrix <- get(file_name)
      
      
      #test_df <-as.data.frame(UNC6_to_UNC15185)
      additive_matrix <- matrix(nrow = nrow(full_matrix), ncol = 8) # empty df with correct row count
      row.names(additive_matrix) <- row.names(full_matrix)
      allele_name <- c("A","B","C","D","E","F","G","H")
      colnames(additive_matrix) <- allele_name
      
      for (allele in allele_name) {
        temp <- full_matrix[, grep(allele, colnames(full_matrix))]
        if (ncol(temp) >= 1) {
          dosage <- 2 * temp[,1]
          if (ncol(temp) > 1) {
            dosage <- dosage + rowSums(temp[, 2:ncol(temp)], na.rm = TRUE)
          }
          additive_matrix[,allele] <- dosage
        }
      }
      
      dir_path1 <- paste0(data_dir,"/additive/",chromosome,"/data/")
      if (!dir.exists(dir_path1)) {
        dir.create(dir_path1, recursive = TRUE)
      }
      
      file_saved = paste0(dir_path1,file_name, ".RData")
      print(file_saved)
      assign(file_name,additive_matrix)
      save(list = file_name, file = file_saved)
      
    }
  }
  
}

create_additve_data("/home/prnpq9/probability_matrices/Genotype",set_chr)
create_additve_data("/home/prnpq9/probability_matrices/updated_Genotype/Genotype",set_chr)
