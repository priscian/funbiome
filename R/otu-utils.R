### Work with CSS data; abbreviate column names and store taxonomy for each column.

#' @export
taxa_names <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

#' @export
taxa_fun <- function(x)
{
  taxa <- sub("^[a-z]__", "", str_split(x, ";")[[1]], perl=TRUE)
  taxa <- gsub("^$", "All", taxa, perl=TRUE)
  length(taxa) <- length(taxa_names)
  names(taxa) <- taxa_names

  as.list(taxa)
}


## Get OTU count before further QC manipulation.
#' @export
GetSmallestTaxonomicGroup <- function(x)
{
  tx <- attr(x, "taxonomy")
  if (is.null(tx))
    return (NULL)

  tail(names(which(!is.na(tx))), 1)
}


#' @export
get_sample_counts_by_taxonomy <- function(m, levels = taxa_names)
{
  summary(factor(unlist(lapply(m, GetSmallestTaxonomicGroup)), levels = levels))
}


## Limit a subset to samples with > 1 observations.
#' @export
create_abundance_data <- function(microbiome, id_variable = "Participant ID")
{
  m <- microbiome
  #m <- dplyr::left_join(m[, .N, by = `Participant ID`][N > 1, 1, with = FALSE], m, by = id_variable)
  m <- dplyr::left_join(m[, .N, by = id_variable][N > 1, 1, with = FALSE], m, by = id_variable)

  m
}


## Remove OTUs whose zero counts >= 98% of all samples.
#' @export
filter_low_abundances <- function(m, zero_ratio_threshold = 0.98)
{
  m <- m[, lapply(.SD,
    function(x) {
      if (!is.null(attr(x, "taxonomy"))) {
        if (quantile(x, zero_ratio_threshold) == 0)
          return (NULL)
      }

      x
    }
  )]

  m
}


## Pick off only the named taxonomic groups from the microbiome samples.
#' @export
GetSamplesByTaxonomicGroupsFactory <- function(taxa) {
  f <- function(x) # Using same name 'taxa' confuses promise evaluation.
  {
    tx <- attr(x, "taxonomy")
    if (!is.null(tx)) {
      if (tail(names(which(!is.na(tx))), 1) %nin% taxa)
        return (NULL)
    }

    x
  }

  f
}

## usage:
# temp <- m[, lapply(.SD, GetSamplesByTaxonomicGroupsFactory(c("class", "phylum")))]


## Log transform the filtered data.
#' @export
transform_abundances <- function(m, trans_fun=function(x) log2(x + 1))
{
  m2 <- data.table::copy(m)
  m2[, 4:ncol(m2) := lapply(m2[, 4:ncol(m2), with=FALSE], trans_fun)] # I.e. use only columns 4+.

  m2
}
