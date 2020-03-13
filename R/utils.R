### Utility functions

#' @export
newdata_from_model <- function(m)
{
  tail(as.data.frame(sapply(rbind(model.frame(m), NA), Hmisc::impute, simplify = FALSE)), 1)
}


## V. http://stackoverflow.com/questions/11121385/repeat-rows-of-a-data-frame
#' @export
rep.data.frame <- function (x, ...) as.data.frame(lapply(x, rep, ...))
