library(dplyr)

#path for the directory of probability Matrix from UNC website
#dir_path <- "/Users/pathumnawarathna/Documents/Research/Simulation/keele/Probability_matrices_from_UNC/test for bug"  
dir_path <- "."

#vector of diploid allele combinations 36 
founder_names <- c("A.A","B.B","C.C","D.D","E.E","F.F","G.G","H.H","A.B","A.C","A.D","A.E","A.F",
                   "A.G","A.H","B.C","B.D","B.E","B.F","B.G","B.H","C.D","C.E","C.F","C.G","C.H","D.E",
                   "D.F","D.G","D.H","E.F","E.G","E.H","F.G","F.H","G.H")

#set of CC lines (number of the CC line)
lines=c(1:63,65,68,70:76,78:83)

#set of chromosomes
set_chr <- c(1:19,"X")

#load Marker file - Marker | chromosome | position.B38
#marker_file <- read.csv(paste(dir_path,"/Markers.csv",sep=""))

#load map file from JAX
#map_file <- read.csv(paste(dir_path,"/mapfile.csv",sep=""))

#load map file from UNC
load(paste0(dir_path,"/snps.megamuga.Rdata"))
map_file <- snps

#change column names of map_file
map_file <- map_file %>% rename("chromosome" = "chr", "position(B38)" = "pos")

#recode chr name 
map_file <- map_file %>%
  mutate(chromosome = recode(chromosome,
                             "chr1" = "1","chr2" = "2","chr3" = "3", "chr4" = "4",
                             "chr5" = "5","chr6" = "6","chr7" = "7", "chr8" = "8",
                             "chr9" = "9","chr10" = "10","chr11" = "11","chr12" = "12", 
                             "chr13" = "13","chr14" = "14","chr15" = "15","chr16" = "16",
                             "chr17" = "17", "chr18" = "18","chr19" = "19","chrP" = "P",
                             "chrX" = "X","chrY" = "Y","chrM" = "M"))


#remove incomplete rows from map_file 
map_file <- map_file[complete.cases(map_file[,c(1:7)]),]

#function to remove replicates and Mit markers (M) in Marker file
remove_rep_and_Mit<-function(df) {
  adj_df <-subset(df,df$chromosome !="M" & !grepl("^rep", df$marker))
  return(adj_df)
}

adj_map_file <- remove_rep_and_Mit(map_file)

# if (!identical(adj_marker_file$marker, adj_map_file$marker)) {
#   stop("Markers are not identical in both content and order")
#   # Reorder adj_map_file to match the order of adj_marker_file$marker
#   adj_map_file <- adj_map_file[match(adj_marker_file$marker, adj_map_file$marker), ]
#   
# }



#remove markers having NA chromosomes
#adj_map_file <- subset(adj_map_file, !is.na(chr))
#adj_marker_file <- subset(adj_marker_file,marker %in% adj_map_file$marker)


#create a vector to save probability matrix and chromosome for each marker
n <- length(adj_map_file$marker)

prob_by_markers <- lapply(seq_len(n), function(i) {
  list(
    matrix_data = matrix(numeric(0), nrow = 0, ncol = 0),
    chr = NA_character_
  )
})

#extact names of makrers from map file
names(prob_by_markers) <- adj_map_file$marker


#assign matrix for each index of vector and assign column and row names for matrices in the vector
for(i in 1:length(adj_map_file$marker)){
  
  prob_by_markers[[i]]$matrix_data <- matrix(nrow = length(lines), ncol = 36)
  colnames(prob_by_markers[[i]]$matrix_data) = founder_names
  rownames(prob_by_markers[[i]]$matrix_data) = rep(NA,length(lines))
  
  prob_by_markers[[i]]$chr <- adj_map_file[i,"chromosome"]
}

# Function to add a string if it doesn't exist already
add_unique_strings <- function(new_strings, string_vector) {
  # Filter only those new_strings that are not in the current vector
  unique_to_add <- new_strings[!(new_strings %in% string_vector)]
  
  string_vector <- c(string_vector, unique_to_add)
  return(string_vector)
}


CC_lines <- c()
for(i in 1:length(lines)){
  index <- lines[i]
  if(index<10){
    pt <- paste0("CC00",index)
  }else{
    pt <- paste0("CC0",index)
  }
  CC_lines[i] <- pt
}

removed_lines <- character(0)
removed_markers <- character(0)

row_count = 1
for (pattern in CC_lines) {
  
  matching_files <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
  print(matching_files)
  
  #generate error if file does not exits
  if(length(matching_files)==0){
    stop(paste("File",pattern,"does not exist !"))
  }
  
  #load probability files 
  raw_CC <- read.csv(matching_files,header = TRUE) 
  
  CC <- remove_rep_and_Mit(raw_CC)
  
  common_markers <- intersect(CC$marker,adj_map_file$marker)
  
  CC <- CC[CC$marker %in% common_markers,]
  
  #check the order of markers of the loaded files match with markers in Marker file 
  # if(!identical(CC$marker,adj_map_file$marker)){
  #   stop(paste("markers",setdiff(adj_map_file$marker,CC$marker),"missing in", pattern, "file"))
  # }
  
  row_sum <- rowSums(CC[,4:ncol(CC)])
  if(all(row_sum==0) | sum(row_sum==0)>10){
    print(paste(pattern,"is removed becouse the file has zero probabilities for all the founder combinations"))
    removed_lines <- add_unique_strings(pattern,removed_lines)
  }else{
    if(any(row_sum==0)){
      removed_markers_temp <- CC$marker[which(row_sum==0)]
      removed_markers <- add_unique_strings(removed_markers_temp,removed_markers)
    }
    
    #check founder allele names same and in same order as input founder names after modifying column names temporarily
    #such that there is a . in between the two founders
    temp_names <- colnames(CC)[4:ncol((CC))]
    expected_cols <- unname(sapply(temp_names,function(x){
      paste0(substr(x,1,1),".",substr(x,2,2))
    }))
    
    
    if (!identical(founder_names, expected_cols)) {
      stop(paste("Founder names do not match the column names in",pattern,"file"))
    }
    
    #assign probabilities to each matrix corresponding to maker by CC lines
    for (j in CC$marker) {
      
      prob_by_markers[[j]]$matrix_data[row_count,] <- as.numeric(CC[CC$marker == j,4:ncol(CC)])
      rownames(prob_by_markers[[j]]$matrix_data)[row_count] <- pattern
      
    }
    row_count = row_count +1
    
  }
  
}
#select the entries in prob_by_markers file corresponding to markers common to both map file and cc probability files 
prob_by_markers <- prob_by_markers[names(prob_by_markers) %in% common_markers]

#remove the set of markers whose probabilities doesn't sum up to 1 from prob_by_markers
prob_by_markers <- prob_by_markers[!names(prob_by_markers) %in% removed_markers]

CC_lines_updated <- CC_lines[!CC_lines %in% removed_lines]

prob_by_markers <- lapply(prob_by_markers,function(x){
  temp_mat <- x$matrix_data
  temp_mat <- temp_mat[c(1:length(CC_lines_updated)),]
  x$matrix_data<- temp_mat
  return(x)
})


#update adj_mapfile to include only the common makers 
adj_map_file <- adj_map_file[adj_map_file$marker %in% common_markers,]

#remove the set of markers whose probabilities doesn't sum up to 1 from adj_map_file
adj_map_file <- adj_map_file[!adj_map_file$marker %in% removed_markers,]

#create list to store bp, chromosome, map, markers, strains, subjects file 
no_chr <- length(set_chr)
sup_files <- lapply(seq_len(no_chr), function(i) {
  list(
    bp = numeric(),
    chromosome = character(),
    map = numeric(),
    markers = character(),
    strains = c("A", "B" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H"),
    subjects = CC_lines_updated
  )
})

save(adj_map_file,file = paste0(dir_path,"/adj_map_file.RData"))

# function to assign values to bp, chromosome, map, markers, strains, subjects for each chromosome and save them 
save_sup_files  <- function(dir_path = dir) {
  names(sup_files)= set_chr
  for (i in set_chr) {
    #assign bp values
    bp <- adj_map_file$`position(B38)`[adj_map_file$chr == i]
    sup_files[[i]]$bp <- bp
    
    #assign chromosome files
    chromosome <- as.character(adj_map_file$chr[adj_map_file$chr == i])
    sup_files[[i]]$chromosome <- chromosome
    
    #assign map files
    map <- adj_map_file$cM[adj_map_file$chr == i]
    sup_files[[i]]$map <- map
    
    #assign maker values
    markers <- adj_map_file$marker[adj_map_file$chr == i]
    sup_files[[i]]$markers <- markers
    
    #assign strains files
    strains <- sup_files[[i]]$strains
    
    #assign subjects values
    subjects <- sup_files[[i]]$subjects
    
    cat(paste0("................Saving files for chr", i,"................"), "\n")
    
    for (folder in c("additive","full","genotype")) {
      
      dir_path1 <- paste0(dir_path, "/Genotype/",folder,"/chr",i,"/")
      
      if (!dir.exists(dir_path1)) {
        dir.create(dir_path1, recursive = TRUE)
      }
      
      
      
      #save bp files
      bp_file_saved = paste0(dir_path1,"bp.RData")
      save(list = "bp", file = bp_file_saved)
      cat(bp_file_saved,"\n")
      
      #save chromosme files
      chr_file_saved = paste0(dir_path1,"chromosome.RData")
      save(list = "chromosome", file = chr_file_saved)
      cat(chr_file_saved,"\n")
      
      #save map files
      map_file_saved = paste0(dir_path1,"map.RData")
      save(list = "map", file = map_file_saved)
      cat(map_file_saved,"\n")
      
      #save marker files
      markers_file_saved = paste0(dir_path1,"markers.RData")
      save(list = "markers", file = markers_file_saved)
      cat(markers_file_saved,"\n")
      
      #save strains files
      strains_file_saved = paste0(dir_path1,"strains.RData")
      save(list = "strains", file = strains_file_saved)
      cat(strains_file_saved,"\n")
      
      #save subjects files
      subjects_file_saved = paste0(dir_path1,"subjects.RData")
      save(list = "subjects", file = subjects_file_saved)
      cat(subjects_file_saved,"\n")
    }
    
    
  }
}


#function to save probability matrices 
save_prob_files <- function(prob_by_markers = prob_list, dir_path = dir){
  for (k in names(prob_by_markers)) {
    assign(k, prob_by_markers[[k]]$matrix_data)
    chr <- prob_by_markers[[k]]$chr
    
    
    chr_dir <- paste0("chr",chr,"/")
    
    dir_path1 <- paste0(dir_path, "/Genotype/full/",chr_dir,"data/")
    
    if (!dir.exists(dir_path1)) {
      dir.create(dir_path1, recursive = TRUE)
    }
    
    
    file_saved = paste0(dir_path1,k, ".RData")
    save(list = k, file = file_saved)
    cat(file_saved,"\n")
    
  }
}


save_sup_files(dir_path)
save_prob_files(prob_by_markers, dir_path)

print("Probability files with default founder order is done")
################################################################################################
#save updated probability matrices with matrix column headers adjusted to fit exactly like keels

#define new order
new_order <- c("A.A","B.B", "C.C" ,"D.D" ,"E.E" ,"F.F", "G.G", "H.H", "B.A", "C.A", "C.B",
               "D.A", "D.B", "D.C", "E.A", "E.B", "E.C", "E.D", "F.A", "F.B", "F.C", "F.D", 
               "F.E", "G.A", "G.B", "G.C", "G.D", "G.E", "G.F", "H.A", "H.B", "H.C", "H.D", "H.E", "H.F", "H.G")

prob_by_markers_updated <- lapply(prob_by_markers,function(x){
  temp_mat <- x$matrix_data
  
  colnames(temp_mat) <- unname(sapply(colnames(temp_mat), function(col_name){
    paste0(substr(col_name,3,3),".",substr(col_name,1,1))
  }))
  
  temp_mat <- temp_mat[,new_order, drop = FALSE]
  
  x$matrix_data<- temp_mat
  return(x)
})

save_sup_files(paste0(dir_path,"/updated_Genotype"))
save_prob_files(prob_by_markers_updated, paste0(dir_path,"/updated_Genotype"))

print("##########################SUMMARY###############################")
print(removed_markers)
print(removed_lines)
