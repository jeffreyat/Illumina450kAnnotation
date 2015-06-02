#'Provides an Illumina450k annotation data access object.
#'Allows for a low memory footprint and efficient data access.
#'
#'@name Illumina450kAnnotationDAO
#'@title Illumina450kAnnotationDAO
#'@field configurer - holds global configuration data in Configurer object
#'@field conn - a connection to a temporary database holding the annotation
#'@field table_name - the name of the table to put data in
#' @export Illumina450kAnnotationDAO
#' @exportClass Illumina450kAnnotationDAO
if(!(require("R6", quietly=T))) {
	install.packages("R6")
	library(R6, quiet=T)
}
Illumina450kAnnotationDAO <- R6Class("Illumina450kAnnotationDAO",
	inherit = AnnotationDAO,
	public = list(
		initialize =  function(configurer) {
			#Create a new DAO object with Illumina450kAnnotationDAO$new(configurer), not by
			#using this method.
			#>configurer -- reference to Configurer
			#<returns nothing
			if(!(require("RSQLite", quietly=T))) {
				install.packages("RSQLite")
				library(RSQLite, quiet=T)
			}
			self$configurer <- configurer
			
			self$table_name <- configurer$getValue('ILLUMINA_ANNOTATION_SQLITE_TABLE')
			
			reg.finalizer(self,
					function(e) {
						if(!is.null(self$conn)){if(dbIsValid(self$conn)){dbDisconnect(self$conn)}}},
					onexit = TRUE)
		}, #end initialize
		
		query_cpg_locations = function(probes) {
			# Return a list of coordinates for CpGs, given a list of probe ids.
			if(!(require("RSQLite", quietly=T))) {
				install.packages("RSQLite")
				library(RSQLite, quiet=T)
			}
			db <- system.file('extdata', 
					self$configurer$getValue('ILLUMINA_ANNOTATION_SQLITE'),
					package='EpiAnalyzer')
			
			self$conn <- dbConnect(SQLite(), dbname=db)
			locs <- dbGetPreparedQuery(self$conn, paste0('select CHR, MAPINFO from ',
							self$table_name,
							' where IlmnID=?'),
					data.frame(probes=probes))
			dbDisconnect(self$conn)
			return(locs)          
		}, # end get_cpg_locations
		
		query_columns = function(columns, constraints=NULL, probe_ids=NULL) {
			# Frequently, the full annotation won't be needed. Just get the
			# columns passed.
			#
			# >columns -- a vector of column names to retrieve
			# >contraints -- such as MAPINFO > 1000000 and CHR = 7
			# >probe_ids -- a vector of probe_ids to keep
			# <returns the contents of the columns asked for
			if(!(require("stringr", quietly=T))) {
				install.packages("stringr")
				library(stringr, quiet=T)
			}
			if(!(require("RSQLite", quietly=T))) {
				install.packages("RSQLite")
				library(RSQLite, quiet=T)
			}
			db <- system.file('extdata', 
					self$configurer$getValue('ILLUMINA_ANNOTATION_SQLITE'),
					package='EpiAnalyzer')
			
			self$conn <- dbConnect(SQLite(), dbname=db)
			sub_table <- NULL
			if(length(constraints)>0) {
				where_clause <- paste0(' where ', constraints)
			} else {
				where_clause <- ''
			}
			
			probe_ids <- data.frame(IlmnID = probe_ids)
			if(length(probe_ids)>0){
				dbGetQuery(self$conn, "attach ':memory:' as mem")
				dbWriteTable(self$conn, "mem.probes", probe_ids)
			}
			
			if(length(probe_ids)>0) {
				if(str_length(where_clause) > 0) {
					probe_clause <- paste0(" and ", self$table_name, 
							".IlmnID = probes.IlmnID")
				} else {
					probe_clause <-  paste0(" where ", self$table_name, 
							".IlmnID = probes.IlmnID")
				}
				
				columns <- paste0(self$table_name, ".", columns)
				col_str <- paste(columns, collapse=", ")
				
				subtable <- dbGetQuery(self$conn, 
						paste0('select ', col_str, 
								' from ', self$table_name,
								', mem.probes ',
								where_clause, probe_clause))
			} else{
				col_str <- paste(columns, collapse=", ")
				subtable <- dbGetQuery(self$conn, paste0('select ', col_str, ' from ', self$table_name, 
								where_clause))
			}
			
			if(length(probe_ids)>0) {
				dbGetQuery(self$conn, 'detach mem')
			}
			
			dbDisconnect(self$conn)
			return(subtable)
		}, #end query_columns
		
		query_probes_by_coordinates = function(coordinates) {
			# Given a set of coordinates, return a data.frame containing info
			# on probes that occure withins those coordinates
			#
			#>coordinates -- a character string in the form of: 
			# chr<num>:<start>-<end>
			#<returns a data.frame containing probe_ids and other
			# information about the probes
			
			if (!(require("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly=T))) {
				biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
				library(TxDb.Hsapiens.UCSC.hg19.knownGene, quiet=T)
			}
			if(!(require("RSQLite", quietly=T))) {
				install.packages("RSQLite")
				library(RSQLite, quiet=T)
			}
			if(!(require("stringr", quietly=T))) {
				install.packages("stringr")
				library(stringr, quiet=T)
			}
			db <- system.file('extdata', 
					self$configurer$getValue('ILLUMINA_ANNOTATION_SQLITE'),
					package='EpiAnalyzer')
			
			self$conn <- dbConnect(SQLite(), dbname=db)
			
			# Get the coordinates of all CpGs in the data.
			loc_matches <- data.frame(str_match_all(coordinates, 
							"chr(\\d+):(\\d+)-(\\d+)"))
			loc_matches$X2 = as.character(paste(loc_matches$X2))
			loc_matches$X3 = as.numeric(paste(loc_matches$X3))
			loc_matches$X4 = as.numeric(paste(loc_matches$X4))
			
			probes <-
					dbGetPreparedQuery(self$conn, 
							paste0("select * from ",
									self$table_name,
									" where CHR = ? ",
									"and MAPINFO between ? and ?"),
							data.frame(chr=loc_matches$X2, 
									loc2=loc_matches$X3, 
									loc2=loc_matches$X4))
			dbDisconnect(self$conn)
			return(probes)
		}, #end query_probes_by_coordinates
		
		query_multi_genes = function () {
			# get probe_ids for probes associated with multiple genes
			
			db <- system.file('extdata', 
					self$configurer$getValue('ILLUMINA_ANNOTATION_SQLITE'),
					package='EpiAnalyzer')
			
			self$conn <- dbConnect(SQLite(), dbname=db)
			# Get the genes associated with each probe.
			gene_col <- self$query_columns(c('IlmnID', 'UCSC_RefGene_Name'))
			
			# Sadly, SQLite cannot split strings, so we have to hack it together.
			# We could come up with in db hack, but memory impact here should be
			# minimal.
			split_genes <- strsplit(as.character(gene_col$UCSC_RefGene_Name), 
					";", fixed=T)
			
			gene_col <- gene_col[sapply(split_genes, function(x) length(unique(x)>1)),]
			
			dbDisconnect(self$conn)
			return(gene_col[,1])
		}, # end query_multi_genes
		
		query_sex_chrom_probes = function() {
			if (!(require("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly=T))) {
				biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
				library(TxDb.Hsapiens.UCSC.hg19.knownGene, quiet=T)
			}
			if(!(require("RSQLite", quietly=T))) {
				install.packages("RSQLite")
				library(RSQLite, quiet=T)
			}
			if(!(require("stringr", quietly=T))) {
				install.packages("stringr")
				library(stringr, quiet=T)
			}
			db <- system.file('extdata', 
					self$configurer$getValue('ILLUMINA_ANNOTATION_SQLITE'),
					package='EpiAnalyzer')
			
			self$conn <- dbConnect(SQLite(), dbname=db)
			
			probes <-
					dbGetQuery(self$conn, 
							paste0("select IlmnID from ",
									self$table_name,
									" where CHR = 'X' or CHR = 'Y' "))
			dbDisconnect(self$conn)
			return(probes)
		} # end query_sex_chrom_probes
	) # end public
) #end Illumina450kAnnotationDAO