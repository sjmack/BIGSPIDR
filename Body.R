#Convert all allele names to 2 fields - Replaces 1 field allele names with NA
#If there are too many alleles that are 1 field, notify user there may be a problem
require(BIGDAWG)


#User chooses a file to be analyzed. This reads the file in as a data frame.
  #Genotype_Data <- read.table(file.choose(), header = TRUE, sep = "\t", quote = "", na.strings = "****", colClasses = "character", check.names = FALSE)

Datafile_Processing <- function(locus, Genotype_Data) {
#Takes every other column and the one after - pairs of 2
  Final_Data <- Genotype_Data[,1:2]
  colnames(Final_Data) <- colnames(Genotype_Data)[1:2] 
#Takes every column pair and runs it though the check function -> Gives a table of the data where the all the alleles are truncated to 2 fields and any 1 field alleles are replaced by NA
  for (x in seq(3,length(Genotype_Data),2)) {
    if (colnames(Genotype_Data[x]) %in% locus) {
      Allele_Columns <- Genotype_Data[,x:(x+1)] ## not a list of lists
      print(paste("Column pairs:", x,(x+1), sep = " "))
      colnames(Allele_Columns) <- colnames(Genotype_Data)[x:(x+1)]
      Final_Data <- cbind(Final_Data, Dataset_Allele_Check_V2(Allele_Columns))
    }
  } 
  Final_Data
}

Dataset_Allele_Check_V2 <- function(Alleles) {
  #Declaring needed variables
  count <- a <- 0
  Temp_List <- apply(Alleles, FUN = GetField, Res = 1, MARGIN = c(1,2))
  Final_Alleles <- data.frame(Alleles, check.names = FALSE) #This will get returned later. We will modify this with the following for loop.
  
  #Takes each column and creates a logical table (T if 1 field, F otherwise) -> Following the logical table, replace data with NA if 1 field, and all other data with 2 field, regardless of initial field count. I.E "12:24" stays "12:24" but "12:52:42" truncates to "12:52"
  for (i in 1:2) {
    comparison <- Alleles[,i] %in% Temp_List[,i]
    count <- sum(length(which(comparison))) + count     #Counts number of 1 field alleles, which show up as TRUE in the comparison table.
    a <- matrix(ifelse(comparison, NA, sapply(Alleles[,i], FUN = GetField, Res = 2)), nrow(Alleles), 1, byrow = FALSE) #a is temporary list for easier replacement of rows.
    Final_Alleles[[i]] <- a
  }
  
  # as.matrix(Final_Alleles)
  
  #Calculates percentage of the data that is 1 field, outputs an integer value denoting how many 1 field alleles were in the data and outputs a percentage.
  percentage <- (count / (nrow(Alleles) * 2))
  print(paste("The number of single field Alleles is:", count, sep = " "))
  print(paste("The percentage of single field Alleles in this column pair is:", percentage, sep = " "))  
  
  #Checks if the percentage of single field alleles is below a certain threshold. This is currently not changable by the user but can be implemented.
  if (percentage > .05) {
    stop("This column pair has too many alleles that are single field.")
  } else {
    print("This column pair is good to go!")
  }
  Final_Alleles
}



  