#combiAnalyzer function script
#by:Livia Tran 
#v 1.3
#4/2/2019

#This script contains code for the function combiAnalyzer(), which recurisvely analyzes 
#BIGDAWG outputs and moves them to the next iteration based on Odds Ratio (OR) improvement.
#OR defines whether a motif is predisposing (OR >1) or protective (OR <1)
#In this case, the threshold for OR improvement is defined at anything greater than 0.1 
#Single variant amino acid positions from the BIGDAWG Output List Object (BOLO) are analyzed
#against a dummy Key Data List Object (KDLO), which is comprised of non-significant values
#i.e p-value = 0.5, OR = 1.0, and CI= 0.5-1.5
#Unassociated motifs are output into an Unassociated Motifs List Object (UMLO) and removed from
#the KDLO
#Variant amino acid positions that move on are combined pairwise; pairs that show OR improvement
#with respect to their single variant amino acids postions move on to form triplets, quadruplets, etc
#until no more improvement is possible 
#Based on the iteration, different pieces of code are called - this is counted by a counter 
#at the beginning of recursion 
#This analyzer currently only supports up to analyzing a maximum of septet combinations
#based on most current protein alignment sequences


#required packages
require(BIGDAWG)
require(gtools)


#loads variantAAtable.rda, which was previously obtained by running variantAA_extractor()
#variantAAtable.rda contains all variant amino acid positions for HLA-A, B, C, DPB1, DRB1, and DQB1
load("variantAAtable.rda")

#begin function script 
combiAnalyzer<-function(loci, myData, KDLO, BOLO, UMLO, counter, motif_list){
  
  #specifies a default motif list if one is not provided 
  if((is.null(motif_list)==TRUE)&(counter==0)){
    motif_list<-c(0,2,3,4,5,6,7)
    print(paste("Motif list has not been provided - function will run until maximal OR is reached"))
  }
      #cat(paste("internal motif_list = ",motif_list,"\n",sep=""))
  
  #BIGDAWG analysis for iteration 0 
  #set output as T for statistical outputs 
  BOLO<-BIGDAWG(myData, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F, Verbose = F)
  
  #unlists all lists in columns in the dataframe 
  BOLO<-data.frame(lapply(as.data.frame(BOLO$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)
  
  #creates dummy_KDLO for comparison to first BOLO ONLY on the 0th iteration 
  if(counter==0){
    #makes dummy KDLO based on previous BOLO 
    dummy_KDLO<-as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F)[rep(seq_len(nrow(as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F))), each=nrow(BOLO)),]
    dummy_KDLO[,1]<-BOLO$Locus
    dummy_KDLO[,2]<-BOLO$Allele
    
    ##MAORI module 
    #finds difference between dummy and BOLO amino acid variants and inputs into new column
    ##dummy comparison only for 0th iteration
    for(i in 1:nrow(BOLO)){
      #finds OR difference between BOLO and dummy ORs -- subs out "-", for a blank, since only evaluating absolute value of OR diff
      #adds difference to new column in BOLO 
      BOLO[i,8]<-gsub("-", "", as.numeric(BOLO[i,]$OR)-as.numeric(subset(subset(dummy_KDLO, grepl(BOLO[i,][[1]], dummy_KDLO[,1])), grepl(BOLO[i,][[2]], subset(dummy_KDLO, grepl(BOLO[i,][[1]], dummy_KDLO[,1]))[,2]))[,3]))[[1]]
    }}
  
  #subsets out binned alleles and any alleles with NA combinations
  if(counter>0){
    BOLO<-subset(BOLO, (BOLO$Allele!="binned") & (!grepl("NA", BOLO$Allele)))}
  
  #MAORI statement for iteration 1
  if(counter==1){
    for(i in 1:nrow(BOLO)){
      BOLO[i,8]<-gsub("-", "", as.numeric(BOLO[i,]$OR)-as.numeric(subset(subset(KDLO, KDLO[,1] %in% strsplit(BOLO[i,][[1]], ":")[[1]][[1]]), subset(KDLO, KDLO[,1] %in% strsplit(BOLO[i,][[1]], ":")[[1]][[1]])$Allele %in% strsplit(BOLO[i,][[2]], "~")[[1]][[1]])$OR))}
  }
  
  #ends function if BOLO is empty 
  if((counter>0) & (nrow(BOLO)==0)){ 
    return(BOLO)}
  
  #MAORI statement for iteration 2+
  #further addition for adding a 9th column for comparison to newly made nth variants to its singular amino acid variant
  if(counter>1){
    for(i in 1:nrow(BOLO)){
      BOLO[i,8]<-gsub("-", "", as.numeric(BOLO[i,]$OR)- as.numeric(subset(subset(KDLO, KDLO[,1] %in% paste(strsplit(BOLO[i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO$Locus, ":")[[1]]))], collapse=":")), subset(KDLO, KDLO[,1] %in% paste(strsplit(BOLO[i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO$Locus, ":")[[1]]))], collapse=":"))$Allele %in%paste(strsplit(BOLO[i,][[2]], "~")[[1]][c(1:length(strsplit(KDLO$Locus, ":")[[1]]))], collapse="~"))$OR))
      BOLO[i,9]<-gsub("-", "", as.numeric(BOLO[i,]$OR)-as.numeric(subset(subset(KDLO_list[[1]], KDLO_list[[1]]$Locus %in% strsplit(BOLO[i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[i,][[1]], ":")))]]), subset(KDLO_list[[1]], KDLO_list[[1]]$Locus %in% strsplit(BOLO[i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[i,][[1]], ":")))]])$Allele %in% strsplit(BOLO[i,][[2]], "~")[[1]][[length(unlist(strsplit(BOLO[i,][[1]], ":")))]])$OR))
    }}

  #subsets out NS values 
  KDLO<-subset(BOLO,BOLO[,7]=="*")

  #statement for returning BOLO if KDLO=0
  if((counter>0) & (nrow(KDLO)==0)){ 
    return(BOLO)}

  #subsets out variants that have not shown >0.1 improvement from their previous variants and
  #singular amino acids 
    if(counter>1){
      
    #subsets out OR differences smaller than 0.1 
    KDLO<-subset(KDLO, KDLO[,9]>0.1)}
  KDLO<-subset(KDLO, KDLO[,8]>0.1)
  
  #statement for returning KDLO if KDLO=0
  if(nrow(KDLO)==0){
    return(KDLO)}
  
  #adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
  KDLO<-unique(rbind(KDLO, subset(BOLO, BOLO$Locus%in%KDLO$Locus)))[mixedorder(row.names(unique(rbind(KDLO, subset(BOLO, BOLO$Locus%in%KDLO$Locus))))),]
  
  #finds unassociated positions from current iteration 
  unassociated_posi<-unique(BOLO$Locus[!BOLO$Locus %in% KDLO$Locus])
  
  #if length(unassociated_posi==0), return KDLO -- this means KDLO and BOLO are the same
  #and max improvement has been reached 
  if(length(unassociated_posi)==0){
    return(KDLO)
  }

  #pair name generation 
  if(counter==0){
    start1<-unique(KDLO$Locus)
    combinames<-sapply(start1, function(x) NULL)
    for(i in 1:(length(start1)-1)){ ## range.x = 1:(N-1)
      for(j in (i+1):length(combinames)){ ## range.y = x+1:N
        if(names(combinames)[[j]]!=start1[[i]]){
          combinames[[i]][[j]]<-paste(start1[[i]],names(combinames)[[j]],sep=":")}}}
    #unlists iter0names and omits NAs to obtain all unique possible pair combinations 
    combinames<-unlist(combinames, use.names = F)[!is.na(unlist(combinames, use.names = F))]
  }
  
  #set start as singular amino acids 
  if(counter>0){
    start1<-unique(unlist(strsplit(KDLO$Locus, ":")))
    combinames<-NULL}
  
  if(counter>0){
    possible_combis<-sapply(unique(KDLO$Locus), function(x) NULL)
    
    #finds possible combinations by pasting names of list with singular amino acids not in that pair 
    for(i in 1:length(possible_combis)){
      possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")}
    
    #splits those triplets up and sorts them numerically to later on eliminate any duplicates 
    for(j in 1:length(unlist(possible_combis))){
      combinames[[j]]<-paste(mixedsort(strsplit(unlist(possible_combis, use.names=F), ":")[[j]], decreasing=F), collapse=":")}
    
    combinames<-unique(mixedsort(combinames))}
  
###subsets combinames by successive unassociated positions
  if(counter==1){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
  }
  
  if(counter==2){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter]])){
     combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
  }
  
  #ends function of no more combination names can be made due to lack of improvement 
  if(length(combinames)==0){
    return(KDLO)
  }

  #triplets
  if(counter==3){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-1]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-1]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-1]][[i]], sep=""), combinames)))}
  }
  if(length(combinames)==0){
    return(KDLO)
  }
  
  #quadruplets
  if(counter==4){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-1]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-1]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-1]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-2]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-2]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-2]][[i]], sep=""), combinames)))}
  }
  
  #ends function of combinames is 0
  if(length(combinames)==0){
    return(KDLO)
  }
  
  #quintets
  if(counter==5){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))} #Vinh's code
    for(i in 1:length(UMLO_list[[counter]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-1]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-1]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-1]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-2]])){
      combinames<-subset(combinames, !grepl(paste("^", UMLO_list[[counter-2]][[i]], sep=""), combinames) & (!grepl(paste(":", UMLO_list[[counter-2]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-3]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-3]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-3]][[i]], sep=""), combinames)))}
  }
  
  #ends function if combinames is 0
  if(length(combinames)==0){
    return(KDLO)
  }
  
  #sextets
  if(counter==6){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-1]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-1]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-1]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-2]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-2]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-2]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-3]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-3]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-3]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-4]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-4]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-4]][[i]], sep=""), combinames)))}
  }
  
  #ends function if combinames is 0
  if(length(combinames)==0){
    return(KDLO)
  }
  
  #septets
  if(counter==7){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-1]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-1]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-1]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-2]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-2]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-2]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-3]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-3]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-3]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-4]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-4]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-4]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-5]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-5]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-5]][[i]], sep=""), combinames)))}
  }
  
  #octets
  if(counter==8){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-1]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-1]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-1]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-2]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-2]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-2]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-3]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-3]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-3]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-4]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-4]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-4]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-5]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-5]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-5]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-6]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-6]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-6]][[i]], sep=""), combinames)))}
  }
  
  #nonetes
  if(counter==9){
    for(i in 1:length(unassociated_posi)){
      combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-1]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-1]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-1]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-2]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-2]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-2]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-3]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-3]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-3]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-4]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-4]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-4]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-5]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-5]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-5]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-6]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-6]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-6]][[i]], sep=""), combinames)))}
    for(i in 1:length(UMLO_list[[counter-7]])){
      combinames<-subset(combinames, (!grepl(paste("^", UMLO_list[[counter-7]][[i]], sep=""), combinames)) & (!grepl(paste(":", UMLO_list[[counter-7]][[i]], sep=""), combinames)))}
    }

  
  #if there are no more combination names after subsetting them based on UMLO_list
  #return function and end
  if(length(combinames)==0){
    return(KDLO)
  }
  
  #df for pairs -- length is number of unique pairs * 2, 
  combidf<-data.frame(variantAAtable[[loci]][,c(1,2)], matrix("", ncol =length(rep(combinames, 2))), stringsAsFactors = F)
  
  #fills in column names 
  colnames(combidf)<-c("SampleID", "Disease", mixedsort(rep(unlist(combinames), 2)))
  
  #observes number of columns for those needed to be pasted together
  cols=c(1:length(strsplit(combinames[[1]], ":")[[1]]))
  
  #[[1]] to contain amino acid combos of TRUE/FALSE
  #[[2]] to contain amino acid combos of FALSE/TRUE
  dfAA<-sapply(1:2, function(x) NULL)
  
  #fills in element names in the lists formed in the above lists 
  for(j in 1:length(dfAA)){
    dfAA[[j]]<-sapply(combinames, function(x) NULL)}
  
  #fills in appropriate position pair combos into dfAA
  for(i in 1:length(combinames)){
    dfAA[[1]][[i]]<-apply(variantAAtable[[loci]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
    dfAA[[2]][[i]]<-apply(variantAAtable[[loci]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  }
  
  #fills into pair_df
  combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
  combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]
  
  #saves each iteration into specified elements in a list in a variable "myData"
  #returns myData
  myData<-list("KDLO"=KDLO, "BOLO"=BOLO, "combidf"=combidf, "UMLO"=unassociated_posi, "combinames"=combinames)
  return(myData)
}


#makes empty lists so results of each iteration may be stored 
BOLO_list<-KDLO_list<-UMLO_list<-list()

#sets motif_list to NULL
motif_list<-NULL

#sets myData, iteration0, to variantAAtable[[loci]]
myData<-variantAAtable[[loci]]

#initiates recursion with stop=FALSE and begins the counter at 0
stop<-FALSE
counter=0

###BEGIN RECURSION -- as long as stop==FALSE, combiAnalyzer will be run until the maximum OR
#is reached, or the end of the motif_list is reached
#the recursive program receives input from combiAnalyzer, where stop=TRUE once the maximum OR 
#is reached, either because the BOLO is empty, the KDLO is empty, or no more combination names
#can be made 
while(stop==FALSE){
  
  #used to inform user what iteration is currently running
  cat(paste(counter,"iteration(s) have been run \n", sep=" "))
  
  #saves each iteration to "interim"
  interim<-combiAnalyzer(loci, myData, BOLO ,KDLO, UMLO, counter, motif_list)
  
  #adds 1 to the counter with each iteration 
  counter=counter+1

  #saves all data to list variables made earlier
  myData<-interim$combidf
  KDLO<-KDLO_list[[counter]]<-interim$KDLO
  BOLO<-BOLO_list[[counter]]<-interim$BOLO
  UMLO<-UMLO_list[[counter]]<-interim$UMLO
  
  #cat(paste("external motif_list = ",motif_list,"\n",sep=""))
  
  if(is.null(nrow(KDLO))==TRUE){
    stop=TRUE
    cat("BIGCAAT -- Maximal OR reached - end of analysis.\n")
  }
  
  if((is.null(nrow(KDLO))==FALSE) & (length(motif_list)!=counter)){
    cat("BIGCAAT -- Dataset is able to be further analyzed - moving on to next iteration.")
  }
  
  if((is.null(nrow(KDLO))==FALSE) & length(motif_list)==counter){
    cat("BIGCAAT -- WARNING: end of motif_list analysis, but further analysis is possible.")
    stop=TRUE
  }
  
  if((is.null(nrow(KDLO))==TRUE) & length(motif_list)==counter){
    cat("BIGCAAT -- End of motif_list analysis - maximal OR has been reached.")}
  
}


