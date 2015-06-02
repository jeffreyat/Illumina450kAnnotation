#'Provides an Illumina450kAnnotation object.
#'This object will hold an annotation in a memory efficient manner.
#'Only the data that is needed will actually be loaded.
#'
#'@name Illumina450kAnnotation
#'@title Illumina450kAnnotation
#'@field configurer - holds global configuration data in Configurer object
#'@field dao - a reference to a DAO for expression data
#' @export Illumina450kAnnotation
#' @exportClass Illumina450kAnnotation
if(!(require("R6", quietly=T))) {
	install.packages("R6")
	library(R6, quiet=T)
}
Illumina450kAnnotation <- R6Class("Illumina450kAnnotation",
	inherit = Annotation,
	public = list(
		configurer = NULL,
		dao = NULL,
		initialize =  function(configurer) {
			#Create a new Illumina450kAnnotation object for holding expression data.
			#But should be done by calling Annotation$new()
			#>annotation -- data.frame with annotation data
			
			self$configurer <- configurer
			
			self$dao <- Illumina450kAnnotationDAO$new(self$configurer)
			
		}, #end initialize
		
		get_columns = function(columns, constraints = NULL, probe_ids = NULL) {
			'This method returns only the columns asked for from the 
					450k annotation.
					
					>columns -- a vector of column names to be returned
					
					>constraints -- optional string containing constraints to place
					on request
					
					>probe_ids -- optional vector of probe_ids to get result for
					
					<returns a data.frame with only the columns asked for of the
					annotation'
			return(self$dao$query_columns(columns, constraints, probe_ids))
		}, #end get_columns
		
		get_cpg_locations = function(probe_ids) {
			'Get CpG locations concordant with given probe ids.
					
					>probe_ids - IlmnID column from 450k annotation file.
					
					<returns a data.frame containing coordinates.
					'
			
			return(self$dao$query_cpg_locations(probe_ids))
		}, # end get_cpg_locations
		
		get_probes_by_coordinates = function(coordinates) {
			' Given a set of coordinates, return a data.frame containing info
			on probes that occure withins those coordinates
			
			>coordinates -- a character string in the form of: 
			 chr<num>:<start>-<end>
			
			<returns a data.frame containing probe_ids and other
			information about the probes'

			return(self$dao$query_probes_by_coordinates(coordinates))
		}, #end get_probes_by_coordinates
		
		get_sex_chrom_probes = function() {
			return(self$dao$query_sex_chrom_probes())
		} # end get_sex_chrom_probes
	) # end public
) #end Illumina450kAnnotation