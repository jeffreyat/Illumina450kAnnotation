#'This class is for getting and setting configuration options. You can 
#'set any text values you like on the fly with 
#'setValue and retrieve them with getValue. They are not saved, but 
#'might be useful during runtime, since the Configurer is available to 
#'all the objects used in this package. 
#' 
#'@name Configurer
#'@title Configurer
#'@field config Reference to a Configurer object.
#' @export Configurer
#' @exportClass Configurer
if(!(require("R6", quietly=T))) {
	install.packages("R6")
	library(R6, quiet=T)
}
Configurer <- R6Class("Configurer",
	public = list(
		config = NULL,
		initialize =  function() {
			self$config <- list(
					DEBUG = 'FALSE',
					REPO = 'http://cran.us.r-project.org',
					EXPRESSION_DATA = '/home/jeff/Downloads/TCGA_BLCA_hMethyl450-2014-08-22/genomicMatrix',
					CGI_ANNOTATION = 'CGI.sqlite',
					REPEAT_ANNOTATION = 'repeats.sqlite',
					RANDOM_SEED = '6433242',
					IDEO_HG19 = 'ideo.RDS',
					ILLUMINA_ANNOTATION = 'HumanMethylation450_15017482_v1-2.csv',
					UCSC_ANNOTATION = 'refGene.txt',
					CR_PROBES = 'crossReactive.txt',
					CORRECT_UCSC_REFGENE_GROUPS = c("TSS1500", "TSS200", "5'UTR", "1stExon"),
					OVERLAP = 'FALSE',
					ILLUMINA_ANNOTATION_RDATA = 'ParsedIlluminaAnnotation.RDS',
					GENE_ANNOTATION_RDATA = 'GeneAnnotation.RDS',
					UCSC_ANNOTATION_FLAGS_RDATA = 'UCSCAnnotationFlags.RDS',
					ILLUMINA2ENTREZ = 'Illumina2Entrez.sqlite',
					ILLUMINA_ANNOTATION_SQLITE = 'Illumina450kAnnotation.sqlite',
					ILLUMINA_ANNOTATION_SQLITE_TABLE = 'Ilmn450kv12',
					THIS_PACKAGE = "RepeatAnnotation"
			) # end list
			if(self$config[['DEBUG']]) {
				self$config[['ILLUMINA_ANNOTATION']] <- 'HumanMethylationShort.csv'
			}
			#callSuper(...)
		}, # end initialize method
		getValue = function(name) {
			'Returns the value of the configuration option asked for using
					a character string.'
			return(self$config[[name]])
		}, #end getValue method
		setValue = function(name, value) {
			'Sets a character string value of the given name in the Configurer.'
			if(name %in% names(config)) {
				self$config[[name]] = value
			} else {
				newlist <- vector(mode='list', length=1)
				names(newlist) = name
				newlist[[name]] = value
				self$config = c(self$config, newlist) 
			}
		} # end setValue method
	) # end public
) # end Configurer