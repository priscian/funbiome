##################################################
### Flexible functions for rapid modeling.
##################################################

## "_M_icrobiome _M_odel"
#' @export
mm <- function(formula, data, cdata=demo, ...)
{
  ## 'formula' is a model formula with a simplified OTU name from the columns of 'data' (e.g. 'colnames(mr2[, 4:NCOL(mr2)])') as the LHS.
  ## 'data' is a data table containing microbiome abundances; here one of 'm[rnt]2' or 'a[rnt]2'.
  ## 'cdata' is a data table or frame containing demographic and clinical model covariates; the default is the data table 'demo'.
  ##   Note that 'pffr()' formulas can't handle variable names with spaces in them.
  ## '...' is any further arguments to be passed to function 'pffr()'.
  ##
  response <- all.vars(formula)[1]
  ## The name of the LHS variable in the formula doesn't matter, but it can't have spaces in it, so fix that:
  f <- eval(substitute(update(formula, `Y` ~ .), list(Y=as.name(make.names(stringr::str_replace_all(response, "(\\s+|\\[|\\]|-)", "_")))))) # "\\u00b7" (∙) doesn't work.
  data.table::setnames(ydata <- dplyr::arrange(data[, c("Participant ID", "ca", response), with=FALSE], `Participant ID`, ca), c(".obs", ".index", ".value"))
  ydata$.obs <- as.numeric(factor(ydata$.obs))

  d <- dplyr::left_join(unique(data[, .(`Participant ID`)]), cdata, by="Participant ID")
  d <- as.data.frame(d); rownames(d) <- NULL
  d$id <- factor(unique(ydata$.obs))
  m <- tryCatch(pffr(f, data=d, ydata=ydata, ...), error=function(e) { cat("ERROR:", conditionMessage(e), "\n"); return (NULL) })

  if (is.null(m)) return (m)

  attr(m, "data_name") <- as.character(as.list(match.call())$data)
  class(m) <- unique(c("MicrobiomeModel", class(m)))

  m
}

## usage:
if (FALSE) {
  m3 <- mm(`s-Enterobacteriaceae;Other.` ~ modB + gabCB + c(sex) + c(race), data=mr2)
  summary(m3)

  ## Use backticks to enclose unsyntactical OTU names (or just always enclose the OTU names with backticks).
  m4 <- mm(`s-[Mogibacteriaceae];All.` ~ modB + gabCB + c(sex) + c(race), data=mr2)
  summary(m4)

  ## Create model fits from a vector of OTU names.
  otus <- colnames(mr2[, 4:NCOL(mr2)])[1:5]
  otus
  f5 <- ~ modB + gabCB + c(sex) + c(race)
  m5 <- sapply(otus, function(x) { f5 <- eval(substitute(update(f5, `Y` ~ .), list(Y=as.name(x)))); mm(f5, data=mr2) }, simplify=FALSE)
  ## N.B. The variable 'm5' contains a list of 'pffr()' models named after their OTU response variables.

  ## Get a list of model fits based on "keystone" OTUs:
  keystonesMap <- otu_taxonomy_map[grep("^(Morax|Coryne|Staphylococ|Haemophil|Streptococ)", names(otu_taxonomy_map), value=TRUE, ignore.case=TRUE)]
  keystoneOtus <- names(keystonesMap)
  m5a <- sapply(intersect(colnames(mr2), keystoneOtus), function(x) { f5 <- eval(substitute(update(f5, `Y` ~ .), list(Y=as.name(x)))); mm(f5, data=mr2) }, simplify=FALSE)
  names(m5a)

  ## IL-8 & CD31
  m6 <- mm(`s-Enterobacteriaceae;Other.` ~ il8B + cd31B + c(sex) + c(race), data=mr2_no_nicu, cdata=demo_il8)
  summary(m6)

  ## Alpha diversity outcome
  m7 <- mm(`PD Whole Tree` ~ modB + gabCB + c(sex) + c(race), data=ar2)
  summary(m7)
}


## Plot function for microbiome models.
#' @export
plot.MicrobiomeModel <- function(x, covariate, covariate_levels, covariate_description="Covariate", subtitle, taxon_level, xlab="Post-Menstrual Age (weeks)", plot=TRUE, ...)
{
  ## 'covariate' is a character string naming a LHS covariate in the formula of 'm'. The effects plot will display the response for all levels of 'covariate'.
  ## 'covariate_levels' gives names to the levels of 'covariate' in the case of 'covariate' being numeric. If missing, the levels are 'names(table(model.frame(x)[[covariate]]))'.
  ## 'covariate_description' is a character string that describes 'covariate' in the plot title.
  ## 'subtitle' is a character string typically naming the swab site of the microbiome abundance samples. If missing, its value is derived from the name of the abundance data set.
  ## 'xlab' is a character string for the x-axis label. The default is "Post-Menstrual Age (weeks)".
  ## 'plot' is a logical value. If TRUE, immediately create the effects plot based on the model 'x' and 'covariate'; if FALSE, return the "ggplot" object.
  ## '...' is further arguments to be passed to function...?
  ##
  newdata <- newdata_from_model(x)
  mfc <- model.frame(x)[[covariate]]
  l <- names(table(mfc))
  if (!is.factor(mfc)) mode(l) <- mode(mfc)

  if (missing(covariate_levels)) {
    covariate_levels <- as.list(l)
    names(covariate_levels) <- l
  }

  response <- all.vars(formula(x))[1]

  if (missing(taxon_level))
    taxon_level <- otu_taxonomy_map[response]

  if (missing(subtitle))
    subtitle <- data_site_map[attr(x, "data_name")]

  ## Expand 'newdata' to have as many rows as no. of levels of the plot covariate.
  newdata <- rep(newdata, length(l))
  newdata[[covariate]] <- l
  newdata$covariate_levels <- factor(l); levels(newdata$covariate_levels) <- covariate_levels

  d <- Reduce(rbind, plyr::alply(newdata, 1,
    function(a)
    {
      yindex.vec <- coef(x)$smterms[[1]]$x
      y <- rep(a, length(yindex.vec))
      y$yindex.vec <- yindex.vec
      p <- predict.gam(x, y, type="response", se.fit=TRUE)

      dataframe(x=yindex.vec, lwr=p$fit - 1.96 * p$se.fit, response=p$fit, upr=p$fit + 1.96 * p$se.fit, group=a$covariate_levels)
    }))

  p <- ggplot() +
    geom_ribbon(data=d, mapping=aes(x=x, ymin=lwr, ymax=upr, fill=group), alpha=0.20) +
    geom_line(data=d, mapping=aes(x=x, y=response, col=group), size=1.0) +
    theme(text=element_text(size=18)) +
    labs(list(title=paste0("Effects Plot of FLM-Fitted Abundance for ", covariate_description), subtitle=subtitle, x=xlab, y=paste0(response, " (", taxon_level, ")"))) +
    labs(caption=paste0("p-value: ", gfp(summary(x)$s.table[paste0(covariate, "(yindex)"), "p-value"]))) + # N.B. Must be separate from the other 'labs()'.
    guides(col=guide_legend(title="Subgroup"), fill=guide_legend(title="Subgroup"))

  if (plot) {
    if (dev.cur() == 1L) # If a graphics device is active, plot there instead of opening a new device.
      dev.new(width=12.5, height=7.3) # New default device of 1200 × 700 px at 96 DPI.

    print(p)
  }
  else
    return (p)
}

## usage:
if (FALSE) {
  m3 <- mm(`s-Enterobacteriaceae;Other.` ~ modB + gabCB + c(sex) + c(race), data=mr2)
  summary(m3)
  p3a <- plot(m3, "gabCB", covariate_levels=list(premature=0, `full-term`=1), covariate_description="Pre- or Full-Term", plot=FALSE)
  p3b <- plot(m3, "modB", covariate_levels=list(vaginal=0, Caesarean=1), covariate_description="Mode of Delivery", plot=FALSE)
  if (dev.cur() == 1L) dev.new(width=12.5, height=7.3)
  print(p3a)
  print(p3b)

  m7 <- mm(`PD Whole Tree` ~ modB + gabCB + c(sex) + c(race), data=ar2)
  plot(m7, "modB", taxon_level="alpha diversity", covariate_levels=list(vaginal=0, Caesarean=1), covariate_description="Mode of Delivery")
  plot(m7, "gabCB", taxon_level="alpha diversity", covariate_levels=list(premature=0, `full-term`=1), covariate_description="Pre- or Full-Term")
}


#' @export
make_spaghetti  <- function(x, ...)
  UseMethod("make_spaghetti")


#' @export
make_spaghetti.pffr <- function(x, ngrid=201, df=20, return_ydata_only=FALSE)
{
  ## plot(x, predict(sm.spline(pd$yindex.vec[1:18], predict(m9, type="lpmatrix")[1:18, ] %*% coef(m9), df=20), x), col="red")

  y <- x$pffr$ydata
  y$id <- factor(y$.obs)
  r <- range(y$.index, na.rm=TRUE)
  xarg <- seq(r[1], r[2], length.out=ngrid)

  response <- predict.gam(x, type="response", se.fit=TRUE)
  ycc <- as.data.frame(y)[as.numeric(rownames(response$fit)), ]
  ycc$row_index <- seq(nrow(ycc))
  ycc <- dataframe(ycc, fit=response$fit, se.fit=response$se.fit)

  if (return_ydata_only)
    return (ycc)

  l <- plyr::dlply(ycc, .(id),
    function(a)
    {
      p <- tryCatch(
        ## { drop(predict(sm.spline(a$yindex.vec, a$fit, df=df), xarg)) }, # Fits spline to fitted data(?)
        ## { xarg <<- a$yindex.vec; a[, 2] }, # Just returns original data.
        ## { xarg <<- a$yindex.vec; a$fit }, # Just returns fitted data.
        { predict.gam(x, newdata=dataframe(id=a$id[1], yindex.vec=xarg), se.fit=TRUE) }, # Returns predictions over entire range of data.
        error = function(e) return (NULL))
      if (is.null(p))
        return (p)

      p <- sapply(p, function(b) { r <- b; is.na(r) <- is.nan(r); r }, simplify=FALSE)

      return (dataframe(id=a$id[1], x=xarg, y=p$fit, lwr=p$fit - 1.96 * p$se.fit, upr=p$fit + 1.96 * p$se.fit, se=p$se.fit))
    })
  d <- Reduce(rbind, l)
  d$id <- factor(d$id)

  attr(d, "OTU") <- split(formula(x))$left_side

  d
}
