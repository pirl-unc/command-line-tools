###### entry point method #########
generate_index_and_run_trace = function( root_directory, entry_point, entry_file, output_file, function_index=NA ){
  
  if ( !is.list(function_index) ) {
    cat("Generating INDEX from nextflow files in: ", root_directory, "\n")
    function_index <- generate_nf_index(root_directory)
    cat("INDEXING completed\n\n")
  } else {
    cat("Using previously created function index.\n\n")
  }
#  return(function_index)
  #print(function_index)
  trace=NULL
  if ( !is.na(entry_point) && !is.na(entry_file) ) {
    cat("Running TRACE from : ", entry_point, " in ", entry_file, "\n\n")
    trace <- trace_nf_path( entry_file, entry_point, function_index=function_index)
    trace %<>% format_trace_for_output()
    if (!is.na(output_file)){
      writeLines(trace, output_file)
      cat("\nTRACE completed, saved to ", output_file, ".\n")
    }
  } else {
    cat("No trace starting point sent. Returning NULL.\n")
  }

  
  return(list(function_index=function_index, trace=trace))
}

#######  supporting methods #########


#######
# trace the path from a given starting point ( workflow or process )
# through each workflow/process called
trace_nf_path <- function( init_file, init_workflow="main", function_index ){
  output_lines <- c()
  this_filename <- strsplit(init_file, "/") %>% {.[[1]][length(.[[1]])]}
  qualified_name <- paste(
    this_filename,
    init_workflow, 
    sep="__" )
  output_lines <- print_this_step( qualified_name, function_index )
  return(output_lines)
}

print_this_step <- function( qualified_name, function_index, level=1, internal_call=NA ){
  cat(rep("  ", (level-1)), 
      (strsplit(qualified_name, "__")[[1]] 
       %>% {paste(.[2], "(", .[1], ")")}), 
         "\n")
#  print(paste("checking ", qualified_name, internal_call))
  if ( !qualified_name %in% names(function_index) ) stop( "Fully qualified workflow/process ", qualified_name, " does not exist in index provided." )
  wfop_object <- function_index[[qualified_name]]
  my_callstack <- paste(wfop_object$file, wfop_object$class, ifelse( is.na(internal_call) | internal_call==wfop_object$name, " ", internal_call ), sep="\t")
  if ( level > 1 ) {
#      print(paste("1", my_callstack, sep=" :: "))
      my_callstack %<>% paste(., paste(rep(" ", level-1), collapse="\t"), sep="\t")
#      print(paste("2", my_callstack, sep=" :: "))
  }
  my_callstack %<>% paste(wfop_object$name, sep="\t")
  
  if ( length(wfop_object$calls) ){
    for( this_subcall in wfop_object$calls ){
      my_callstack %<>% c( print_this_step( this_subcall$qualified_name, function_index, level+1, this_subcall$internal_name ) )
    }
  }
  return(my_callstack)
}

format_trace_for_output <- function( trace ){
  #determine max number of columns in trace
  num_columns <- mapply(function(x){
    return(length(strsplit(x, "\t")[[1]]))
  }, trace) %>% max()
  #add appropriate empty cells to end of each line
  trace %<>% mapply( function( line ){
    line %<>% {strsplit(., split="\t")[[1]]} %>%
      {c(.,rep(" ", num_columns - length(.)))} %>%
      paste(collapse="\t")
    return(line)
  }, .)
  #add column headers
  trace %<>% c(paste("FILE", "ACTION", "INTERNAL_NAME", "BASE_METHOD", paste("LEVEL", 1:(num_columns-4), collapse="\t", sep="_"), sep="\t"), .)
  return(trace)
}

a <- function( v ){
  index_output <<- c(index_output, v)
  cat(v)
}
####### PARSING .nf FILES for INCLUDES #########

###########  
# check all .nf files in the directory sent and all subdirectories
# create index of all workflows/processes with callstacks for 
# workflows and processes invoked by each workflow
generate_nf_index <- function( my_directory ){
  
  nf_files <- dir(my_directory, pattern="*\\.nf", recursive=T, full.names=TRUE)# %>% .[grep(".nf$", .)]
#  print(nf_files)
  if ( length(nf_files) == 0 ) 
    stop("The directory sent doesn't have any .nf files to process.")
  
  #confirm no duplicate filenames
  names_only <- nf_files %>% gsub("^[a-zA-Z_0-9]+/", "", .)
  if ( length(names_only) != length(unique(names_only)) )
    stop("The following filenames are duplicated ... sorry, this script won't work properly::\n", nf_files[duplicated(names_only)])
  
  #generate index of all nf functions ( processes and workflows )
  index_output <<- c()
  function_list <- list()
  for ( nf_file in nf_files ){
#    cat(paste("Indexing", nf_file), "\n")
    split_path <- strsplit(nf_file, "/")[[1]]
    processing_filename <- split_path[ length(split_path) ]
    
    a(paste0("\nProcessing FILE :: ", processing_filename, "\n"))
    
    namespace <- ifelse( length(split_path) > 1, split_path[length(split_path)-1], NA )
    fl <- file(nf_file)
    fl_lines <- readLines( fl )
    close(fl)
    
    all_includes <- list()
    my_includes <- list()
    in_script_block <- FALSE
    for ( line in fl_lines ) {
      line <- trimws(line, which="left")
      if ( line == "" | grepl( "^#", line) ) next
      #skip script blocks
      if ( grepl("^[\"\']{1,3}[\t ]*$", line) ){
        in_script_block %<>% !.
        next
      } else if ( in_script_block ) {
        next
      }
      
      # check for include statment
      this_include <- check_for_include(line)
      #  this is an include block line
      if ( length( this_include ) > 1 || this_include == TRUE ) {
        if ( is.logical(this_include) ) {
          #multi-line include block ... get next line
          next
        }
        if ( length(this_include$modules) > 0 ) {
          #there are new modules to be included ...
          for ( m in this_include$modules ) {
            my_includes[ m[ "internal_name" ] ] <- m[ "actual_name" ]
          }
        }
        if ( !is.na(this_include[ "filename" ]) ) {
          #we've got the full list of modules in this include block
          #push them onto the all_includes list
          for ( inc in names(my_includes) ) {
            
            a(paste0("INCLUDE ::", paste( my_includes[ inc ], this_include[ "filename" ], sep=" from " ),"\n"))
            
            all_includes[ inc ] <- paste( this_include[ "filename" ], my_includes[ inc ], sep="__" )
          }
          my_includes <- list()
        }
        next
      }
      # check for process or workflow
      process <- check_for_process_or_workflow(line, processing_filename)
      if ( is.list(process) ) {
        qualified_process <- paste( process$file, process$name, sep="__" )
        
        a(paste("Found PROCESS ::", process$name, "\n"))
        
        if( qualified_process %in% names(function_list) ){
          warning("PROCESS/WORKFLOW :: ", processing_filename, "::", process, " exists in another namespace ( from ", function_list[[qualified_process]]$file, " )")
        } else {
          function_list[[ qualified_process ]] <- process
        }
        next
      }
      
      # check for call to process or workflow, push onto call stack for the current workflow/process
      call <- check_for_call( line, all_includes, processing_filename )
      if ( is.list( call ) ) {
        
        a(paste("\tcalls ::", call[[1]], "\n"))

        function_list[[ length(function_list) ]]$calls %<>% append(list(call))
        next
      }
      
    }
#    break
  }
  return(function_list)
}
#now


# include { 
#   rscript;  
#   generic_process_1_out as coxph_os_ranks;
# } from './generic/generic.nf'
#  value: character vector as [1] internal name used, [2] include-filename__actual-name

check_for_include = function( line ){
  module_inclusion = NULL
  from_file = NULL
  # 1- include {
  # 2- include { something } from somewhere.nf
  # 3- include { something as something_else } from somewhere.nf
  if ( grepl( "^include[ ]*\\{", line ) ) {
    #check if there is a module block on this inclusion line
    if ( grepl("\\{[\t ]*[a-zA-Z0-9_-]+", line) ) {
      #strip leading characters / white space
      module_inclusion <- gsub("[^\\{]+\\{[\t ]*", "", line)
    } else {
      #this is an empty include statement ... modules to follow
      return( TRUE )
    }
  # 4- something;
  # 5- something as something_else;
  } else if ( grepl( "^[\t ]*[a-zA-Z0-9_-]+([ ]+as)|([ ]*;)", line) ) {
    #strip off leading spaces or tabs from module statement
    module_inclusion <- gsub( "^[\t ]*", "", line )
  }
  #strip following content from module_inclusion found
  if ( !is.null(module_inclusion) ) {
    # strip following ;, } ... 
    module_inclusion %<>% gsub("[\t ]*\\}.*$", "", .) %>% gsub("[ ]*;$", "", .)
  }
  # 6- ... } from somewhere.nf
  if ( grepl( "\\}[ ]*from", line ) ){
    from_file <- gsub(".*\\}[ ]*from[ ]+", "", line)  %>% gsub("[ ]+$", "", .)
  }
  rtn <- list( modules=list(), filename=NA )
  if ( !is.null( module_inclusion ) ) {
    modules <- strsplit(module_inclusion, "[\t ]*;[\t ]*")[[1]]
    for ( m in modules ) {
      name_qualifier <- strsplit(m, "[ ]+(AS|as)[ ]+")[[1]]
      actual_name <- name_qualifier[1]
      internal_name <- ifelse(length(name_qualifier) > 1, name_qualifier[2], actual_name)
      rtn$modules[[ length(rtn$modules) +1 ]] <- c( internal_name = internal_name, actual_name = actual_name )
    }
  }
  
  if ( !is.null( from_file ) ) {
    include_filename <- strsplit(from_file, "/") %>% {.[[1]][length(.[[1]])]} %>% gsub("('|\")", "", .)
    rtn[ "filename" ] = include_filename
  }
  
  if ( !is.na( rtn$filename ) || length(rtn$modules) > 0 ){
    return( rtn )
  } else {
    return( FALSE )
  }
}

check_for_process_or_workflow = function( line, processing_filename ){
  process <- grep("^(process|workflow) [a-zA-Z0-9_]*[ ]?\\{", line, value=TRUE)
  if ( length(process) ) {
    f_class <- ifelse( grepl("^process", process), "process", "workflow" )
    process %<>% gsub("^(process|workflow)[ ]+", "", .) %>% 
      gsub("[ ]*\\{.*", "", .)
    if( process == "" ) process = "main"
    return(list(class=f_class, name=process, file=processing_filename, calls=list()))
  }
  return(FALSE)
}

check_for_call = function( line, includes, processing_filename ){
  call <- grep("^[A-Za-z0-9_]+[ ]*\\(", line, value=TRUE)
  if ( length(call) ) {
    call_name <- gsub( "[ ]*\\(.*$", "", call )
    #there are a few calls that we might catch which aren't actually calls - exclude them here ...
    if( call_name %in% c("path", "val", "tuple", "if", "evaluate") ) return(FALSE)
    # find this call in includes, or if not there, append current filename since we know it must be in this namespace
    qualified_name <- ifelse( call_name %in% names(includes), includes[[ call_name ]], paste( processing_filename, call_name, sep="__") )
    return( list( internal_name=call_name, qualified_name=qualified_name ) )
  }
  return(FALSE)
}


########### functional code - at the end so it can be run from command line and have all methods pre-created #########

if ("magrittr" %in% rownames(installed.packages()) == FALSE) {
  install.packages('magrittr', repos = "http://cran.us.r-project.org")
}
library(magrittr)

####### init variables #########

root_directory <- file.path(getwd() ) #, "viral")
out_directory <- root_directory
entry_file <- NA #"main.nf"
entry_point <- NA #"main" #use main to enter in the "default" or unnamed workflow in entry_file, otherwise use name of workflow
nf_index <- NA
run_stats <- FALSE

####### override defaults with command line arguments if any were sent ########
args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  root_directory <- args[1]
  out_directory <- args[2]
  #check directory exists
  if ( !dir.exists(root_directory) ) {
    stop( "The directory sent does not appear to exist." )
  }
  if (length(args) >= 4) {
    entry_file <- args[3]
    entry_point <- args[4]
  }
}

output_filename <- paste("trace",paste0(entry_point, ".tsv"), sep="_" )

#index_rds_path <- file.path(root_directory, "nf_index.RDS")

#if ( file.exists( index_rds_path ) ) nf_index <- readRDS(index_rds_path)


# still finding a problem parsing some of the files from the nkc project ... 


index_output <- c()
fi <- generate_index_and_run_trace( root_directory, entry_point, entry_file, file.path(out_directory,output_filename), function_index = nf_index )

### save index output to trace annotation file ...


# optionally save index to disk
#saveRDS(fi$function_index, index_rds_path)

# optionally run stats on the trace results
if ( run_stats ) {
  trace_df <- read.table(file.path(root_directory, output_filename), sep="\t", header = TRUE)
  length(unique(trace_df$LEVEL_1))
  length(trace_df$LEVEL_1[trace_df$LEVEL_1 != " "])
}

