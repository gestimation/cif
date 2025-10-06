#' @title Visualize time-to-event outcomes and intercurrent events
#' @description
#' Estimate and plot survival curves using the Kaplan–Meier estimator or
#' cumulative incidence curves under competing risks using the Aalen–Johansen estimator.
#' Returns a \code{survfit}-compatible object and, by default, draws a publication-ready plot via \pkg{ggsurvfit}.
#'
#' @details
#' \strong{Estimation:}
#' \itemize{
#'   \item \code{outcome.type = "SURVIVAL"}: Kaplan–Meier estimator with Greenwood-type variance.
#'   \item \code{outcome.type = "COMPETING-RISK"}: Aalen–Johansen estimator for CIF of \code{code.event1}
#'         using IPCW for the censoring distribution. The returned \code{surv} corresponds to \code{1 - CIF}.
#' }
#' \strong{Confidence intervals:}
#' Constructed on the probability scale with the specified \code{conf.type}.
#' If \code{conf.type \%in\% c("none","n")}, the plot suppresses CI bands.
#'
#' \strong{Plotting:}
#' By default, the function calls an internal helper \code{call_ggsurvfit()} which adds
#' confidence bands, risk table, censoring marks, and optional intercurrent-event marks.
#' For CIF display, set \code{ggsurvfit.type = "risk"}.
#'
#' @param formula A model formula specifying the outcome and (optionally) \code{strata()}.
#' @param data A data frame containing variables in \code{formula}.
#' @param weights Optional name of the weight variable in \code{data}. Weights must be nonnegative; strictly positive is recommended.
#' @param subset.condition Optional character expression to subset \code{data} before analysis.
#' @param na.action Function to handle missing values (default: \code{\link[stats]{na.omit}}).
#' @param code.event1 Integer code of the event of interest (default \code{1}).
#' @param code.event2 Integer code of the competing event (default \code{2}).
#' @param code.censoring Integer code of censoring (default \code{0}).
#' @param outcome.type \code{"SURVIVAL"} (KM) or \code{"COMPETING-RISK"} (AJ).
#' @param conf.int Two-sided confidence level (default \code{0.95}).
#' @param error Character specifying variance type used internally.
#'   For \code{"SURVIVAL"} typically \code{"greenwood"}; for \code{"COMPETING-RISK"} pass options supported by \code{calculateAalenDeltaSE()} (e.g., \code{"aalen"}, \code{"delta"}, \code{"none"}).
#' @param conf.type Transformation for CI on the probability scale (default \code{"arcsine-square root"}).
#' @param report.survfit.std.err If \code{TRUE}, report SE on the log-survival scale (survfit's convention). Otherwise SE is on the probability scale.
#' @param report.ggsurvfit If \code{TRUE} (default), draw a \pkg{ggsurvfit} plot.
#' @param ggsurvfit.type \code{NULL} (survival) or \code{"risk"} (display \code{1 - survival} i.e. CIF).
#' @param addConfidenceInterval Logical add \code{add_confidence_interval()} to plot (default \code{TRUE}).
#' @param addRiskTable Logical add \code{add_risktable(risktable_stats="n.risk")} to plot (default \code{TRUE}).
#' @param addCensorMark Logical add \code{add_censor_mark()} to plot (default \code{TRUE}).
#' @param shape.censor.mark Integer point shape for censor marks (default \code{3}).
#' @param size.censor.mark Numeric point size for censor marks (default \code{2}).
#' @param addIntercurrentEventMark Logical; overlay user-specified time marks per strata (default \code{TRUE}).
#' @param intercurrent.event.time Named list of numeric vectors (names must match or be mapped to strata labels).
#' @param shape.intercurrent.event.mark Integer point shape for intercurrent-event marks (default \code{16}).
#' @param size.intercurrent.event.mark Numeric point size for intercurrent-event marks (default \code{2}).
#' @param label.strata Optional character vector of labels for strata.
#' @param label.x,label.y Axis labels (defaults: \code{"Time"}, \code{"Survival probability"}).
#'   If \code{ggsurvfit.type="risk"} and \code{label.y} is unchanged, it is internally set to \code{"1 - survival probability"}.
#' @param lims.x,lims.y Numeric length-2 vectors for axis limits (defaults: \code{NULL}, \code{c(0,1)}).
#' @param font.family,font.size Plot theme controls (defaults: \code{"sans"}, \code{14}).
#' @param legend.position Legend position: \code{"top"}, \code{"right"}, \code{"bottom"}, \code{"left"}, or \code{"none"} (default \code{"top"}).
#'
#' @returns A \code{survfit} object. For \code{outcome.type="SURVIVAL"}, \code{$surv} is the survival function.
#' For \code{outcome.type="COMPETING-RISK"}, \code{$surv} equals \code{1 - CIF} for \code{code.event1}.
#' Standard error and CIs are provided per \code{conf.type}. Note that some methods for \code{survfit} (e.g., \code{residuals.survfit}) may not be supported.
#'
#' @seealso \code{\link{cif_reg}} for log-odds product models of CIFs; \pkg{ggsurvfit} for plotting helpers.
#'
#' @examples
#' library(cif)
#' data(diabetes.complications)
#' survfit_by_group <- cif_curve(Event(t,epsilon) ~ fruitq, data = diabetes.complications,
#'                     outcome.type='COMPETING-RISK', error='delta', ggsurvfit.type = 'risk',
#'                     label.y = 'CIF of diabetic retinopathy', label.x = 'Years from registration')
#'
#' @importFrom ggsurvfit ggsurvfit add_confidence_interval add_risktable add_censor_mark
#' @importFrom ggplot2 theme_classic theme element_text labs lims geom_point aes
#' @importFrom Rcpp sourceCpp
#' @useDynLib cif, .registration = TRUE
#' @export
cif_curve <- function(formula,
                           data,
                           weights = NULL,
                           subset.condition = NULL,
                           na.action = na.omit,
                           code.event1 = 1,
                           code.event2 = 2,
                           code.censoring = 0,
                           outcome.type = "SURVIVAL",
                           conf.int = 0.95,
                           error = "greenwood",
                           conf.type = "arcsine-square root",
                           report.survfit.std.err = FALSE,
                           report.ggsurvfit = TRUE,
                           ggsurvfit.type = NULL,
                           addConfidenceInterval = TRUE,
                           addRiskTable = TRUE,
                           addCensorMark = TRUE,
                           addIntercurrentEventMark = TRUE,
                           label.x = "Time",
                           label.y = "Survival probability",
                           label.strata = NULL,
                           lims.x = NULL,
                           lims.y = c(0, 1),
                           shape.censor.mark = 3,
                           size.censor.mark = 2,
                           intercurrent.event.time = NULL,
                          #intercurrent.event.time = list("fruitq1=0" = c(1, 3), "fruitq1=1" = c(2, 6)),
                           shape.intercurrent.event.mark = 16,
                           size.intercurrent.event.mark = 2,
                           font.family = "sans",
                           font.size = 14,
                           legend.position = "top") {
  outcome.type <- check_outcome.type(outcome.type)
  out_readSurv <- readSurv(formula, data, weights, code.event1, code.event2, code.censoring, subset.condition, na.action)
  error <- check_error(error, outcome.type)
  check_label.strata(out_readSurv, label.strata)

  if (outcome.type == "SURVIVAL") {
    out_km <- calculateKM(out_readSurv$t, out_readSurv$d, out_readSurv$w, as.integer(out_readSurv$strata), error)
    out_km$std.err <- out_km$surv * out_km$std.err
    out_ci <- calculateCI(out_km, conf.int, conf.type, conf.lower=NULL)

    if (!all(as.integer(out_readSurv$strata) == 1) & is.null(label.strata)) {
      names(out_km$strata) <- levels(as.factor(out_readSurv$strata))
    } else if (!all(as.integer(out_readSurv$strata) == 1)) {
      names(out_km$strata) <- label.strata
    }
    if (report.survfit.std.err) {
      out_km$std.err <- out_km$std.err / out_km$surv
    }

    survfit_object <- list(
      time = out_km$time,
      surv = out_km$surv,
      n = out_km$n,
      n.risk = out_km$n.risk,
      n.event = out_km$n.event,
      n.censor = out_km$n.censor,
      std.err = out_km$std.err,
      upper = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$upper,
      lower = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$lower,
      conf.type = conf.type,
      call = match.call(),
      type = "kaplan-meier",
      method = "Kaplan-Meier"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) {
      survfit_object$strata <- out_km$strata
    }
    class(survfit_object) <- c("survfit")
  } else {
    out_aj <- calculateAJ(out_readSurv)
    out_aj <- readStrata(out_readSurv, out_aj, label.strata)
    if (any(as.integer(out_readSurv$strata) != 1)) {
      n <- table(as.integer(out_readSurv$strata))
      rep_list <- mapply(rep, n, out_aj$strata1, SIMPLIFY = FALSE)
      n.risk <- do.call(c, rep_list) - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    } else {
      n <- length(out_readSurv$strata)
      n.risk <- n - out_aj$n.cum.censor - out_aj$n.cum.event1 - out_aj$n.cum.event2
    }
    out_aj$std.err <- calculateAalenDeltaSE(out_aj$time1, out_aj$aj1, out_aj$n.event1, out_aj$n.event2, n.risk, out_aj$time0, out_aj$km0, out_aj$strata1, error)
    out_aj$surv <- 1 - out_aj$aj1
    out_ci <- calculateCI(out_aj, conf.int, conf.type, conf.lower=NULL)
    if (report.survfit.std.err) {
      out_aj$std.err <- out_aj$std.err / out_aj$surv
    }

    survfit_object <- list(
      time = out_aj$time1,
      surv = out_aj$surv,
      n = n,
      n.risk = n.risk,
      n.event = out_aj$n.event1,
      n.censor = out_aj$n.censor,
      std.err = out_aj$std.err,
      upper = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$upper,
      lower = if (is.null(conf.type) | conf.type == "none" | conf.type == "n") NULL else out_ci$lower,
      conf.type = conf.type,
      call = match.call(),
      type = "Aalen-Johansen",
      method = "Aalen-Johansen"
    )
    if (any(as.integer(out_readSurv$strata) != 1)) {
      survfit_object$strata <- out_aj$strata1
    }
    class(survfit_object) <- c("survfit")
  }

  if (report.ggsurvfit) {
    out_ggsurvfit <- call_ggsurvfit(
      survfit_object = survfit_object,
      out_readSurv   = out_readSurv,
      ggsurvfit.type = ggsurvfit.type,
      conf.type      = conf.type,
      addConfidenceInterval = addConfidenceInterval,
      addRiskTable   = addRiskTable,
      addCensorMark  = addCensorMark,
      shape.censor.mark = shape.censor.mark,
      size.censor.mark  = size.censor.mark,
      addIntercurrentEventMark = addIntercurrentEventMark,
      intercurrent.event.time  = intercurrent.event.time,
      shape.intercurrent.event.mark = shape.intercurrent.event.mark,
      size.intercurrent.event.mark  = size.intercurrent.event.mark,
      label.x = label.x, label.y = label.y,
      lims.x = lims.x, lims.y = lims.y,
      font.family = font.family, font.size = font.size,
      legend.position = legend.position
    )
    print(out_ggsurvfit)
  }
  return(survfit_object)
}


#' Plot survival or cumulative incidence curves with ggsurvfit
#'
#' @param survfit_object A \code{survfit} object.
#' @param out_readSurv (optional) List returned by your \code{readSurv()} to auto-set x limits.
#' @param ggsurvfit.type Character; NULL (survival) or "risk" for cumulative incidence display.
#' @param conf.type Character; same as used when constructing CI (e.g., "none", "n", "arcsine-square root").
#' @param addConfidenceInterval Logical; add CI via \code{add_confidence_interval()}.
#' @param addRiskTable Logical; add risk table via \code{add_risktable(risktable_stats="n.risk")}.
#' @param addCensorMark Logical; add censor marks via \code{add_censor_mark()}.
#' @param shape.censor.mark Integer; shape for censor marks. Defaults to 3.
#' @param size.censor.mark Numeric; size for censor marks. Defaults to 2.
#' @param addIntercurrentEventMark Logical; add user-specified marks at strata-specific times.
#' @param intercurrent.event.time Named list of numeric vectors.
#' @param shape.intercurrent.event.mark Integer; fixed shape for intercurrent marks across strata. Defaults 16.
#' @param size.intercurrent.event.mark Numeric; size for intercurrent marks. Defaults 2.
#' @param label.x,label.y Axis labels. If \code{label.y = "Survival probability"} (default) and
#'   \code{ggsurvfit.type == "risk"}, it is automatically replaced with \code{"1 - survival probability"}.
#' @param lims.x Numeric length-2; x limits. If NULL and \code{out_readSurv} given, uses \code{c(0,max(out_readSurv$t))}.
#' @param lims.y Numeric length-2; y limits.
#' @param font.family,font.size Theme controls.
#' @param legend.position "top","right","bottom","left" or "none".
#' @return A \code{ggplot} object.
call_ggsurvfit <- function(
    survfit_object,
    out_readSurv = NULL,
    ggsurvfit.type = NULL,
    conf.type = "arcsine-square root",
    addConfidenceInterval = TRUE,
    addRiskTable = TRUE,
    addCensorMark = TRUE,
    shape.censor.mark = 3,
    size.censor.mark = 2,
    addIntercurrentEventMark = TRUE,
    intercurrent.event.time = list(),
    shape.intercurrent.event.mark = 16,
    size.intercurrent.event.mark = 2,
    label.x = "Time",
    label.y = "Survival probability",
    lims.x = NULL,
    lims.y = c(0, 1),
    font.family = "sans",
    font.size = 14,
    legend.position = "top"
){
  if (is.null(lims.x)) {
    if (!is.null(out_readSurv) && !is.null(out_readSurv$t)) {
      lims.x <- c(0, max(out_readSurv$t))
    }
  }

  label.y.corrected <- if (identical(ggsurvfit.type, "risk") && identical(label.y, "Survival probability")) {
    "1 - survival probability"
  } else {
    label.y
  }

  check_ggsurvfit(
    survfit_object = survfit_object,
    lims.x = lims.x, lims.y = lims.y,
    ggsurvfit.type = ggsurvfit.type,
    addConfidenceInterval = addConfidenceInterval,
    addCensorMark = addCensorMark,
    addIntercurrentEventMark = addIntercurrentEventMark,
    shape.censor.mark = shape.censor.mark,
    shape.intercurrent.event.mark = shape.intercurrent.event.mark
  )

  select_y_axis <- function(sfobj) {
    if (identical(ggsurvfit.type, "risk")) ggsurvfit(sfobj, type = "risk") else ggsurvfit(sfobj)
  }

  if (conf.type %in% c("none","n") || length(survfit_object$strata) > 2) {
    survfit_object_ <- survfit_object
    survfit_object_$lower <- survfit_object_$surv
    survfit_object_$upper <- survfit_object_$surv

    p <- select_y_axis(survfit_object_) +
      theme_classic() +
      theme(
        legend.position = legend.position,
        axis.title = element_text(size = (font.size + 2), family = font.family),
        axis.text  = element_text(size = font.size, family = font.family),
        legend.text = element_text(size = font.size, family = font.family)
      ) +
      labs(x = label.x, y = label.y.corrected) +
      lims(x = lims.x, y = lims.y) +
      theme(legend.position = "top")

    if (isTRUE(addConfidenceInterval)) p <- p + add_confidence_interval()
    if (isTRUE(addRiskTable))         p <- p + add_risktable(risktable_stats = c("n.risk"))
    if (isTRUE(addCensorMark))        p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)

  } else {
    if (length(survfit_object$time) != length(survfit_object$lower))
      stop("time and upper/lower required for ggsurvfit are different lengths")
    if (length(survfit_object$time) != length(survfit_object$n.risk))
      stop("time and n.risk used required ggsurvfit are different lengths")

    p <- select_y_axis(survfit_object) +
      theme_classic() +
      theme(
        legend.position = legend.position,
        axis.title = element_text(size = (font.size + 2), family = font.family),
        axis.text  = element_text(size = font.size, family = font.family),
        legend.text = element_text(size = font.size, family = font.family)
      ) +
      labs(x = label.x, y = label.y.corrected) +
      lims(x = lims.x, y = lims.y)

    if (isTRUE(addConfidenceInterval)) p <- p + add_confidence_interval()
    if (isTRUE(addRiskTable))         p <- p + add_risktable(risktable_stats = c("n.risk"))
    if (isTRUE(addCensorMark))        p <- p + add_censor_mark(shape = shape.censor.mark, size = size.censor.mark)
  }

  if (isTRUE(addIntercurrentEventMark) && length(intercurrent.event.time)) {

    align_marks_keys <- function(fit, marks) {
      canon_str <- function(x) sub("^.*=", "", as.character(x))
      if (is.null(marks) || length(marks) == 0) return(marks)
      fit_names <- if (is.null(fit$strata)) "(all)" else names(fit$strata)
      fit_key   <- canon_str(fit_names)
      out <- list()
      for (k in names(marks)) {
        k_can <- canon_str(k)
        j <- match(k_can, fit_key)
        if (is.na(j)) next
        out[[ fit_names[j] ]] <- marks[[k]]
      }
      out
    }

    make_mark_df <- function(fit, marks, extend = TRUE) {
      if (is.null(marks) || length(marks) == 0) return(NULL)
      marks2 <- align_marks_keys(fit, marks)
      strata_names <- if (is.null(fit$strata)) "(all)" else names(fit$strata)

      out <- lapply(names(marks2), function(st_lab) {
        tt <- marks2[[st_lab]]
        if (length(tt) == 0) return(NULL)
        fit_st <- if (is.null(fit$strata)) {
          if (!identical(st_lab, "(all)")) return(NULL)
          fit
        } else {
          i <- match(st_lab, strata_names); if (is.na(i)) return(NULL)
          fit[i]
        }
        sm <- summary(fit_st, times = tt, extend = extend)
        y  <- if (identical(ggsurvfit.type, "risk")) 1 - sm$surv else sm$surv
        data.frame(strata = st_lab, time = sm$time, y = y, stringsAsFactors = FALSE)
      })
      do.call(rbind, out)
    }

    mark_df <- make_mark_df(survfit_object, intercurrent.event.time)
    if (!is.null(mark_df) && nrow(mark_df) > 0) {
      p <- p + geom_point(
        data = mark_df,
        aes(x = time, y = y, group = strata, colour = strata),
        inherit.aes = FALSE,
        shape = shape.intercurrent.event.mark,
        size  = size.intercurrent.event.mark,
        show.legend = FALSE
      )
    }
  }
  return(p)
}

calculateAJ <- function(data) {
  out_km0 <- calculateKM(data$t, data$d0, data$w, as.integer(data$strata), "none")
  km0 <- get_surv(data$t, out_km0$surv, out_km0$time, as.integer(data$strata), out_km0$strata)
  ip.weight <- (data$d0 == 0) * ifelse(km0 > 0, 1 / km0, 0)
  d1_ipw <- as.matrix(data$w * data$d1 * ip.weight)

  aj1 <- time1 <- integer(0)
  n.cum.event1 <- n.cum.event2 <- n.cum.censor <- numeric(0)
  n.event1 <- n.event2 <- n.censor <- numeric(0)
  strata1 <- integer(0)
  strata_vec <- as.integer(data$strata)

  for (level in sort(unique(strata_vec))) {
    idx <- which(strata_vec == level)

    sub_t  <- data$t[idx]
    sub_d0 <- data$d0[idx]
    sub_d1 <- data$d1[idx]
    sub_d2 <- data$d2[idx]
    sub_d1_ipw <- d1_ipw[idx, , drop = FALSE]

    o <- order(sub_t)
    sub_t  <- sub_t[o]
    sub_d0 <- sub_d0[o]
    sub_d1 <- sub_d1[o]
    sub_d2 <- sub_d2[o]
    sub_d1_ipw <- sub_d1_ipw[o, , drop = FALSE]

    not_atrisk <- outer(sub_t, sub_t, FUN = ">=")
    sub_aj1 <- as.vector(not_atrisk %*% sub_d1_ipw) / length(sub_t)
    sub_n.censor <- as.vector(not_atrisk %*% as.matrix(sub_d0))
    sub_n.event1 <- as.vector(not_atrisk %*% as.matrix(sub_d1))
    sub_n.event2 <- as.vector(not_atrisk %*% as.matrix(sub_d2))

    keep <- !duplicated(rev(sub_t))
    keep <- rev(keep)

    u_t  <- sub_t[keep]
    u_aj1 <- sub_aj1[keep]
    u_nc  <- sub_n.censor[keep]
    u_ne1 <- sub_n.event1[keep]
    u_ne2 <- sub_n.event2[keep]

    oo <- order(u_t)
    u_t  <- u_t[oo]
    u_aj1 <- u_aj1[oo]
    u_nc  <- u_nc[oo]
    u_ne1 <- u_ne1[oo]
    u_ne2 <- u_ne2[oo]

    inc_nc  <- c(u_nc[1],  diff(u_nc))
    inc_ne1 <- c(u_ne1[1], diff(u_ne1))
    inc_ne2 <- c(u_ne2[1], diff(u_ne2))

    time1 <- c(time1, u_t)
    aj1   <- c(aj1, u_aj1)

    n.cum.censor <- c(n.cum.censor, u_nc)
    n.cum.event1 <- c(n.cum.event1, u_ne1)
    n.cum.event2 <- c(n.cum.event2, u_ne2)

    n.censor <- c(n.censor, inc_nc)
    n.event1 <- c(n.event1, inc_ne1)
    n.event2 <- c(n.event2, inc_ne2)

    strata1 <- c(strata1, length(u_t))
  }

  list(
    time1 = time1,
    aj1 = aj1,
    n.event1 = n.event1,
    n.event2 = n.event2,
    n.censor = n.censor,
    n.cum.event1 = n.cum.event1,
    n.cum.event2 = n.cum.event2,
    n.cum.censor = n.cum.censor,
    strata1 = strata1,
    time0 = out_km0$time,
    km0 = out_km0$surv
  )
}

get_surv <- function(predicted.time, estimated.surv, estimated.time, predicted.strata=NULL, estimated.strata=NULL) {
  predicted.surv <- numeric(length(predicted.time))
  strata_start <- c(1, head(cumsum(estimated.strata), -1) + 1)
  strata_end <- cumsum(estimated.strata)
  if (any(is.na(predicted.time)))
    stop("Invalid predicted time variable. NA values included")

  for (i in seq_along(predicted.time)) {
    t <- predicted.time[i]
    strata_size <- estimated.strata[predicted.strata[i]]

    if (is.null(estimated.strata)|all(is.na(estimated.strata))) {
      #      time_until_t <- estimated.time[estimated.time <= t]
      time_until_t <- estimated.time[estimated.time < t]
      if (length(time_until_t) > 0) {
        time_index <- which.max(time_until_t)
        predicted.surv[i] <- estimated.surv[time_index]
      } else {
        predicted.surv[i] <- 1
      }
    } else if (strata_size > 0|all(is.na(estimated.strata))) {
      strata_indices <- strata_start[predicted.strata[i]]:strata_end[predicted.strata[i]]
      strata_time <- estimated.time[strata_indices]
      strata_surv <- estimated.surv[strata_indices]

      #      time_until_t <- strata_time[strata_time <= t]
      time_until_t <- strata_time[strata_time < t]

      if (length(time_until_t) > 0) {
        time_index <- which.max(time_until_t)
        predicted.surv[i] <- strata_surv[time_index]
      } else {
        predicted.surv[i] <- 1
      }
    } else {
      predicted.surv[i] <- NA
    }
  }
  return(predicted.surv)
}

