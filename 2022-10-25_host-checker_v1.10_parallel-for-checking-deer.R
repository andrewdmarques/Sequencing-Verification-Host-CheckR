# Load libraries
library(RColorBrewer)
library(kableExtra)
library(foreach)    # To run distance calculations in parallel
library(doParallel) # To run distance calculations in parallel
library(ggplot2)

# User defined variables. __________________________________________________________________________
make_database <- "no"
sample_size   <- 10000
prefix <- "2022-10-25_weiss-and-28s-pilot_100k"
sample_file <- "2022-10-25_weiss-and-28s-host-check-samples.csv"

# File architecture prep.  __________________________________________________________________________
if(!file.exists("Host")){system("mkdir Host")}
if(!file.exists("Sample")){system("mkdir Sample")}
if(!file.exists("Result")){system("mkdir Result")}


# Function prep.  __________________________________________________________________________
unix <- function(command){
  print(command)
  # Write a bash script that will run the blast.
  filepath01 <- '.'
  p <- c('#!/bin/bash', command)
  writeLines(p, file.path(filepath01, 'temp.script'))
  system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
  comm <- paste0(file.path(filepath01, 'temp.script'))
  # Run the bash script.
  system(comm)
}


# Database Prep  __________________________________________________________________________

# Define the samples that should be ine database.
host_accession  <- c("GCF_000001405.39", "GCF_002102435.1", "GCF_018350175.1", "NC_045512.2")
host_organism     <- c("Human", "Deer", "Cat", "SARS-CoV-2")

# Generate the reference file for host samples.
col <- c("host_organism", "host_accession", "file_location", "command_unzip", "command_combine", "command_make_db")
ref_host <- data.frame(matrix(NA, nrow = length(host_accession), ncol = length(col)))
colnames(ref_host)      <- col
ref_host$host_accession <- host_accession
ref_host$host_organism    <- host_organism
ref_host$file_location  <- paste0("Host/", ref_host$host_organism, "/")
ref_host$command_unzip  <- paste0("unzip ", ref_host$file_location, ref_host$host_accession, ".zip")
ref_host$command_combine <- paste0("cat ", ref_host$file_location, "Database/*.fna > ", ref_host$file_location, "Database/", ref_host$host_organism, ".fasta")
ref_host$command_make_db <- paste0("makeblastdb -in ", ref_host$file_location, "Database/", ref_host$host_organism, ".fasta -dbtype nucl -parse_seqids")

# Create the database. 
if(make_database == "yes"){
  
  # Execute the commands for making the database.
  for(i in 1:length(ref_host$host_organism)){
    if(ref_host$host_organism[i] == "SARS-CoV-2"){
      # Make a database direcotry and put the fasta files there.
      if(!file.exists(paste0(ref_host$file_location[i], "Database"))){system(paste0("mkdir ", ref_host$file_location[i], "Database"))}
      unix(paste0("cp ", ref_host$file_location[i], ref_host$host_accession[i], ".fasta ", ref_host$file_location[i], "Database/"))
      unix(paste0("mv ", ref_host$file_location[i], "Database/", ref_host$host_accession[i], ".fasta ", ref_host$file_location[i], "Database/", ref_host$host_organism[i], ".fasta "))
    } else {
      # Unzip.
      unix(ref_host$command_unzip[i])
      # Move the files to the correct location.
      unix(paste0("mv ncbi_dataset ./", ref_host$file_location[i]))
      unix(paste0("mv README.md ./", ref_host$file_location[i]))
      # Make a database direcotry and put the fasta files there.
      if(!file.exists(paste0(ref_host$file_location[i], "Database"))){system(paste0("mkdir ", ref_host$file_location[i], "Database"))}
      unix(paste0("cp ", ref_host$file_location[i], "ncbi_dataset/data/", ref_host$host_accession[i], "/*.fna ", ref_host$file_location[i], "Database/"))
      # Combine all the fasta files into one fasta file.
      unix(ref_host$command_combine[i])
      # Remove all of the original fasta files.
      unix(paste0("rm -r ", ref_host$file_location[i], "Database/*.fna"))
    }
    # Make the database.
    unix(ref_host$command_make_db[i])
  }
}


# Sample Prep  __________________________________________________________________________

# Determine the samples that should be checked.
sam <- read.csv(sample_file)
for(i in 1:length(sam$VSP)){
  if(sam$VSP[i] == "all"){
    temp <- list.files(paste0("/data/SARS-CoV-2/sequencing/", sam$run[i]))
    temp <- subset(temp, grepl("R1", temp)) 
    col <- c("VSP", "run")
    sam2 <- data.frame(matrix(NA, nrow = length(temp), ncol = length(col)))
    colnames(sam2)      <- col 
    sam2$VSP <- temp
    sam2$run <- sam$run[i]
    sam <- sam2
    sam <- sam[- grep("ndetermined", sam$VSP),]
    sam <- sam[- grep("VSP9", sam$VSP),]
    sam$VSP <- gsub("-.*","",sam$VSP)
    # Comment this out, it is for troubleshooting.
    # sam <- tail(sam, n=30)
    # sam <- head(sam, n=15)
    
  }
}

# Generate the reference file for the samples.
col <- c("VSP", "run", "original_location", "file_location", "error", "sequences_checked")
for(i in 1:length(ref_host$host_organism)){
  temp <- paste0(tolower(ref_host$host_organism[i]), "_hits")
  col <- c(col, temp)
}
ref_sample <- data.frame(matrix(NA, nrow = length(sam$VSP), ncol = length(col)))
colnames(ref_sample) <- col
ref_sample$VSP <- sam$VSP
ref_sample$run <- sam$run
ref_sample$error <- "None"
ref_sample$sequences_checked <- sample_size
ref_sample$original_location <- paste0("/data/SARS-CoV-2/sequencing/", ref_sample$run, "/", ref_sample$VSP, "*R1*")
ref_sample$file_location <- paste0("Sample/", ref_sample$VSP, "/")

# Iterate through each sample of interest and prepare it for analysis.

main_function <- function(ref_sample,i,prefix,ref_host){
  # Make the directory for each sample.
  if(!file.exists(paste0(ref_sample$file_location[i]))){system(paste0("mkdir ", ref_sample$file_location[i]))}
  # Copy the fastq file to this location.
  system(paste0("cp ", ref_sample$original_location[i], " ./", ref_sample$file_location[i]))
  # Convert the fastq.gz file to a fasta file.
  system(paste0("seqtk seq -a ", ref_sample$file_location[i], "*R1*.fastq.gz > ", ref_sample$file_location[i], "sample.fasta"))
  # Take the first sequences (specified in sample_size variable) and place them in a new file.
  system(paste0("sed -n 1,", as.character(2*sample_size), "p ", ref_sample$file_location[i], "sample.fasta > ", ref_sample$file_location[i], "/sample_", as.character(sample_size), ".fasta"))
  # Blast the sample against each of the host databases created.
  for(j in 1:length(ref_host$host_organism)){
    # Blast.
    system(paste0("blastn -db ", ref_host$file_location[j], "Database/", ref_host$host_organism[j], ".fasta -query ", ref_sample$file_location[i], "sample_", as.character(sample_size), ".fasta -out ", ref_sample$file_location[i], "sample-", as.character(sample_size), "_", ref_host$host_organism[j], "-db_result.out"))
    # # Determine the number of sequences that DID NOT match the host blasted against. ###NOTE: These commented lines will make the counts but they will have duplicates of the reads.
    # unix(paste0("grep -o '***** No hits found *****' ", ref_sample$file_location[i], "sample-", as.character(sample_size), "_", ref_host$host_organism[j], "-db_result.out > Sample/temp.csv"))
    # # Record the number of sequences that DID match the host blasted against.
    # temp <- data.frame(readLines("Sample/temp.csv"))
    # colnames(temp) <- "temp"
    # ref_sample[ref_sample$VSP==ref_sample$VSP[i], paste0(tolower(ref_host$host_organism[j]), "_hits")] <- sample_size - length(temp$temp)
  }
  
  # Determine the closest match for each sequence (if a read scores more closely to Human than Deer, then only count it toward Human)
  # Determine all the reads that were blasted.
  system(paste0("grep Query= ", ref_sample$file_location[i], "sample-", as.character(sample_size), "_Human-db_result.out > ", ref_sample$file_location[i], "temp_sample-", as.character(sample_size), "_blasted.temp"))
  fasta_name <- read.csv(paste0(ref_sample$file_location[i], "temp_sample-", as.character(sample_size), "_blasted.temp"), header = FALSE)
  # Make a blank data frame that will have all scores for most similar hits.
  col2 <- c("name", "closest_host")
  for(j in 1:length(ref_host$host_organism)){
    temp <- paste0(tolower(ref_host$host_organism[j]), "_hits")
    col2 <- c(col2, temp)
  }
  fasta_name2 <- data.frame(matrix(NA, nrow = length(fasta_name$V1), ncol = length(col2)))
  colnames(fasta_name2) <- col2
  fasta_name2$name <- fasta_name$V1
  fasta_name2$name <- gsub("Query= ","",fasta_name2$name)
  fasta_name2$name <- gsub(" .*","",fasta_name2$name)
  # Populate the data frame with the score for the closest hits of each read for each host checked.
  for(j in 1:length(ref_host$host_organism)){
    system(paste0("grep -B 5 -A 2 significant ", ref_sample$file_location[i], "sample-", as.character(sample_size), "_", ref_host$host_organism[j], "-db_result.out > ", ref_sample$file_location[i], "temp_sample-", as.character(sample_size), "_", ref_host$host_organism[j], "-db_result.temp"))
    # Only perform the following function if file contains data. No data means that no reads mapped to that host.
    # Make a file that lists all the fasta names from those that had hits.
    system(paste0("grep Query= ", ref_sample$file_location[i], "temp_sample-", as.character(sample_size), "_", ref_host$host_organism[j], "-db_result.temp", " > ", ref_sample$file_location[i], "temp_sample-hits.temp"))
    info = file.info(paste0(ref_sample$file_location[i], "temp_sample-hits.temp"))
    if(info$size > 0){
      temp <- read.csv(paste0(ref_sample$file_location[i], "temp_sample-hits.temp"), header = F)
      temp$V1 <- gsub("Query= ","",temp$V1)
      temp$V1 <- gsub(" .*","",temp$V1)
      for(k in 1:length(temp$V1)){
        for(l in 1:length(fasta_name2$name)){
          if(fasta_name2$name[l] == temp$V1[k]){
            fasta_name2[fasta_name2$name==fasta_name2$name[l], paste0(tolower(ref_host$host_organism[j]), "_hits")] <- "match"
          }
        }
      }
    }
  }
  fasta_name2[is.na(fasta_name2)] <- ""
  for(j in 1:length(fasta_name2$name)){
    for(k in length(col2):1){
      if(fasta_name2[j,k] == "match"){
        fasta_name2[j,2] <- col2[k]
      }
    }
  }
  # Record the counts to the ref_sample dataframe.
  temp <- data.frame(table(fasta_name2$closest_host))
  temp$Var1 <- as.character(temp$Var1)
  temp$Freq <- as.numeric(temp$Freq)
  for(j in 1:length(temp$Var1)){
    if(temp$Var1[j] != ""){
      ref_sample[ref_sample$VSP==ref_sample$VSP[i], temp$Var1[j]] <- temp$Freq[j]
    }
  }
  # Clean up the environment to take up less space.
  system(paste0("rm -r ", ref_sample$file_location[i], "sample.fasta"))
  system(paste0("rm -r ", ref_sample$file_location[i], "temp_*"))
  # Save an intermediate copy of the sample reference file.
  ref_sample[is.na(ref_sample)] <- 0
  # write.csv(ref_sample, paste0("Result/", prefix, "_reference_sample.csv"), row.names = T)
  return(ref_sample[i,])
}



# Prepare for Parallel.
n.cores <- parallel::detectCores() - 1
n.cores <- 22
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK",outfile="")
doParallel::registerDoParallel(cl = my.cluster)

# Determine the nearest neighbor using a parallel function.
ref_sample2 <- data.frame(matrix(NA, nrow = 0, ncol = length(colnames(ref_sample))))
colnames(ref_sample2) <- colnames(ref_sample)
ref_sample1 <- ref_sample
ref_sample2 <- foreach(l = 1:length(ref_sample$VSP),.combine = rbind) %dopar% {
  main_function(ref_sample1,l,prefix,ref_host)
}
ref_sample <- ref_sample2
write.csv(ref_sample, paste0("Result/", prefix, "_reference_sample.csv"), row.names = T)

# Graph the data  __________________________________________________________________________
data_row <- c(ref_host$host_organism, "Other")
dat1 <- data.frame(matrix(NA, nrow = length(data_row), ncol = length(ref_sample$VSP)))
colnames(dat1) <- ref_sample$VSP
rownames(dat1) <- data_row

# Populate the data frame.
for(i in 1:length(ref_sample$VSP)){
  temp <- ref_sample[i,]
  temp <- temp[ , -which(names(temp) %in% c("VSP","run", "original_location", "file_location", "error", "sequences_checked"))]
  temp <- lapply(temp,as.numeric)
  temp <- data.frame(temp)
  temp <- cbind(temp, other_hits = sample_size - rowSums(temp))
  temp <- data.frame(t(temp))
  temp$description <- row.names(temp)
  rownames(temp) <- NULL
  colnames(temp) <- c("data", "description")
  for(j in 1:length(data_row)){
    for(k in 1:length(temp$description)){
      if(temp$description[k] == gsub("-", ".", paste0(tolower(data_row[j]), "_hits"))){
      dat1[j,i] <- temp$data[k]
    }}
  }
}

# Determine the "other" category.
dat2 <- data.frame(t(dat1))
for(i in 1:length(dat2$Human)){
  if(dat2$Other[i] < 0){
    if("SARS-CoV-2" %in% data_row){
      dat2$SARS.CoV.2[i] <- dat2$SARS.CoV.2[i] + dat2$Other[i]
    }
    dat2$Other[i] <- 0
  }
}
dat3 <- dat2[order(dat2$Human),]    # Do this to plot sorted by number of human reads.
# dat3 <- dat2[order(dat2$Deer),]    # Do this to plot sorted by number of deer reads.
dat3 <- data.frame(t(dat3))

# create color palette:
color <- brewer.pal(length(data_row), "Spectral")

# Transform this data in Proportion
dat3 <- apply(dat3, 2, function(x){x/sample_size})

# Make a stacked barplot--> it will be in %!
pdf(file = paste0("Result/", prefix, "_summary-graph.pdf"),   # The directory you want to save the file in
    width = 6.5*1.618,
    height = 6.5)

par(mar = c(5, 4, 4, 8) + 0.1)
barplot(dat3,                     # Draw barplot with properly aligned legend
        col = color,
        border = "white",
        legend.text = TRUE, 
        args.legend = list(x = "topright",
                           inset = c(- 0.15, 0)),
        ylab="Proportion of Reads",
        las = 2, cex.names = 0.25)
dev.off()

# Output a nice table.
colnames(dat2) <- data_row
dat2 <- apply(dat2, 2, function(x){x/sample_size})
dat2 <- apply(dat2, 2, function(x){round(x,2)})
kbl(dat2)

table <- dat2 %>%
  kbl(caption = "Proportion of Reads") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  kable_styling(latex_options="scale_down")

save_kable(table, paste0("Result/", prefix, "_summary-table.pdf"))

# Save the csv file.
dat4 <- data.frame(t(dat3))
dat4$VSP <- rownames(dat4)
dat4 <- dat4[order(dat4$VSP),]
dat4 <- dat4[ , c("VSP", names(dat4)[names(dat4) != "VSP"])]
colnames(dat4) <- gsub("SARS.CoV.2", "SARS-CoV-2", colnames(dat4))
write.csv(dat4, paste0("Result/", prefix, "_summary.csv"), row.names = F)

# For the deer paper, compare to human VSP summarized. ################################################################################
dat4 <-read.csv(paste0("Result/", prefix, "_summary.csv"))
rownames(dat4) <- dat4$VSP 
dat4 <- dat4[ , -which(names(dat4) %in% c("VSP"))]

# Condense all of the non-deer samples into one row.
special_sample <- c("VSP3514","VSP3515","VSP3516","VSP3517","VSP3573","VSP3574","VSP3575")
dat5 <- dat4
dat5$VSP <- rownames(dat5)
dat5 <- dat5[order(dat5$VSP),] 
dat5 <- dat5[ , -which(names(dat5) %in% c("VSP"))]
dat5[nrow(dat5) + 1,] = 0
dat5$filter <- F
rownames(dat5)[length(dat5$Human)]<-"Human Average"
# Iterate through all of the human sampes and add them to the human average.
for(i in 1:length(dat5$Human)){
  if(rownames(dat5)[i] %in% special_sample){
    print(rownames(dat5)[i])
    dat5$filter[i] <- T
  }else{
    if(i != length(dat5$Human)){
    dat5$Human[length(dat5$Human)] <- dat5$Human[length(dat5$Human)] + dat5$Human[i]
    dat5$Deer[length(dat5$Human)] <- dat5$Deer[length(dat5$Human)] + dat5$Deer[i]
    dat5$Cat[length(dat5$Human)] <- dat5$Cat[length(dat5$Human)] + dat5$Cat[i]
    dat5$SARS.CoV.2[length(dat5$Human)] <- dat5$SARS.CoV.2[length(dat5$Human)] + dat5$SARS.CoV.2[i]
    dat5$Other[length(dat5$Human)] <- dat5$Other[length(dat5$Human)] + dat5$Other[i]
    }else{
      dat5$filter[i] <- T
    }
  }
}
# Average the average.
denominator <- (length(dat5$Human)-length(special_sample)-1)
dat5$Human[length(dat5$Human)] <- dat5$Human[length(dat5$Human)]/denominator
dat5$Deer[length(dat5$Deer)] <- dat5$Deer[length(dat5$Deer)]/denominator
dat5$Cat[length(dat5$Cat)] <- dat5$Cat[length(dat5$Cat)]/denominator
dat5$SARS.CoV.2[length(dat5$SARS.CoV.2)] <- dat5$SARS.CoV.2[length(dat5$SARS.CoV.2)]/denominator
dat5$Other[length(dat5$Other)] <- dat5$Other[length(dat5$Other)]/denominator

# Remove all non special sample rows.
dat5 <- subset(dat5, dat5$filter == T)
dat5 <- dat5[ , -which(names(dat5) %in% c("filter"))]
dat6 <- dat5

# dat3 <- data.frame(t(dat4))
dat5 <- data.frame(t(dat5))
dat5 <- apply(dat5, 2, function(x){x})

# create color palette:
color <- brewer.pal(length(data_row), "Spectral")

pdf(file = paste0("Result/", prefix, "_summary-graph_average-human.pdf"),   # The directory you want to save the file in
    width = 6.5*1.618,
    height = 6.5)
par(mar = c(7,5,3,9))
barplot(dat5,                     # Draw barplot with properly aligned legend
        col = color,
        border = "white",
        legend.text = TRUE, 
        args.legend = list(x = "topright",
                           inset = c(-0.2, 0)),
        ylab="Proportion of Reads",
        las = 2, 
        cex.names = 0.9)
dev.off()


linearize <- function(data_frame){
  r <- rownames(data_frame)
  c <- colnames(data_frame)
  
  # Make a blank data frame.
  col <- c('dim1','dim2', 'value')
  linear <- data.frame(matrix(NA, nrow = length(r)*length(c), ncol = length(col)))
  colnames(linear ) <- col
  
  x <- 1
  y <- 1
  for(i in 1:length(linear$dim1)){
    linear$dim1[i] <- r[x]
    linear$dim2[i] <- c[y]
    linear$value[i] <- data_frame[x,y]
    y <- y + 1
    if(y > length(c)){
      y <- 1
      x <- x + 1
      
    }
  }
  return(linear)
}

# Make a grouped bar plot ########################################################
dat6 <- dat6[ , which(names(dat6) %in% c("Human", "Deer"))]

# Lineagize the data.
dat7 <- linearize(dat6)
colnames(dat7) <- c("Sample", "Read Origin", "Proportion of Reads")

# Plot the data
plot_grouped <- ggplot(dat7, aes(fill=`Read Origin`, y=`Proportion of Reads`, x=Sample)) + 
  geom_bar(position='dodge', stat='identity') +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.4)) +
  scale_fill_brewer(palette="Dark2") + 
  theme(axis.title.x = element_blank()) 

plot_grouped

ggsave(paste0("Result/", prefix, "_host-read-origin.pdf"), plot_grouped)
