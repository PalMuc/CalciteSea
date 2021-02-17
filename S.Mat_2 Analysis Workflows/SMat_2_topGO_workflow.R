> library("topGO")


					#### IMPOR DATA AND CREATE topGO OBJECT ####


#import all sequences IDs and associated GO terms

> database_GO <- readMappings(file = "PATH TO GO_TERMS FILE") 
> gene_database <- names(geneID2GO)


#import differentially expressed genes IDs (one per line)

> genes_of_interest<- read.table("interestinggenes.txt",header=FALSE)
> genes_of_interest <- as.character(genesOfInterest$V1


#find genes of interest in the gene_database

> genes_list <- factor(as.integer(gene_database %in% genes_of_interest))
> names(genes_list) <- gene_database


#create topGO object for analysis. Onthologies: BP (biological process), CC (cellular component) and MF (molecular function).

> GO_data <- new("topGOdata", description="Project_name", ontology="BP", allGenes=genes_list,  annot = annFUN.gene2GO, gene2GO = gene_database)


#check that topGO object was successfully created

> myGOdata 



					#### PERFROM ENRICHMENT ANALYSIS ####
	
#run Fisher exact test	
					
> result_Fisher_Test <- runTest(GO_data, algorithm="weight01", statistic="fisher")


# obtain significant results. Number of most significant X terms to be displayed can be changed by changing the topNodes argument in the function below.

> allRes <- GenTable(GO_data, classicFisher = resultFisher, orderBy = "result_Fisher_Test", ranksOf = "classicFisher", topNodes = 10)


# Visualize position of enriched GO Terms in the GO hierarchy.

> showSigOfNodes(GO_data, score(result_Fisher_Test), firstSigNodes = 5, useInfo ='all')


> printGraph(GO_data, result_Fisher_Test, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)


#Obtain genes (from the database file) associated with significantly enriched GO Terms.


> GO_terms_of_interest = c("")

> genes_associated_to_GO_term <- genesInTerm(GO_data, GO_terms_of_interest)

> for (i in 1:length(GO_terms_of_interest))

   {
       term <- GO_terms_of_interest[i]
       gene_with_term_term <- enes_associated_to_GO_term[term][[1]]
       gene_with_term_term <- paste(gene_with_term_term, collapse=',')
       print(paste("Term",term,"genes:",gene_with_term_term))
     }