#' Functions to calculate boundary nutrient concentrations
#' @name boundaries
#' @return a named vector of estimated nutrient boundaries on the scale of \code{xvar}
#' @examples
NULL

#' @describeIn boundaries A wrapper function to estimate boundaries using all methods
#' @param data Appropriate water quality data for estimating nutrient status boundaries
#' @param yvar The name of the response variable, typically EQR. Defaults to y.u
#' @param xvar The name of the pressure variable, e.g. total phosphorus. Defaults to x.u
#' @param status_boundaries A named vector of status boundaries on the scale of \code{yvar}. 
#'      Default contains 0.8, 0.6 and 0.4 for the High-Good, Good-Moderate and Moderate-Poor boundaries 
#'
#' @export
boundaries_all <- function(data,
                       yvar = "y.u",
                       xvar = "x.u",
                       status_boundaries = c("HG" = 0.8, "GM" = 0.6, "MP" = 0.4),
                       class_var = "status"){
  
  OLS <- boundaries_ols(data, yvar, xvar, status_boundaries)
  RMA <- boundaries_RMA(data, yvar, xvar, status_boundaries)
  MM <- boundaries_MM(data, yvar, xvar, status_boundaries)
  AdjQ <- boundaries_AdjQ(data, yvar, xvar, class_var)
  Median <- boundaries_medians(data, yvar, xvar, class_var)
  q75 <- boundaries_q75(data, yvar, xvar, class_var)
  
  out <- data.frame(rbind(OLS, RMA, MM, AdjQ, Median, q75))
  boundary_names <- colnames(out)
  out$Method = rownames(out)
  rownames(out) <- NULL
  out <- out[, c(ncol(out), 1:ncol(out)-1)]
  return(out)
  #tidyr::gather_(out, "Boundary", "Value", gather_cols = boundary_names)
  
}

#' @describeIn boundaries Estimates nutrient boundaries using ordinary-least-squares
#' @inheritParams boundaries_lm
#' @export
boundaries_ols <- function(data,
                          yvar = "y.u",
                          xvar = "x.u",
                          status_boundaries = c("HG" = 0.8, "GM" = 0.6, "MP" = 0.4)) {
  
  f <- as.formula(paste0(yvar, "~", xvar))
  lm1 <- lm(f, data = data)
  cfs <- coef(lm1)
  
  boundaries <- (status_boundaries - cfs[1]) / cfs[2]
  
  return(boundaries)
  
}

#' @describeIn boundaries Estimates nutrient boundaries using Ranged Major Axis regression.
#' @inheritParams boundaries_lm
#' @export
boundaries_RMA <- function(data,
                             yvar = "y.u",
                             xvar = "x.u",
                             status_boundaries = c("HG" = 0.8, "GM" = 0.6, "MP" = 0.4)) {
  
  f <- as.formula(paste0(yvar, "~", xvar))
  lmod <-
    lmodel2::lmodel2(f,
                     range.y = "interval",
                     range.x = "interval",
                     nperm = 99,
                     data = data)
  
  
  RMA.int <- lmod$regression.results[[2]][[4]]
  RMA.slope <- lmod$regression.results[[3]][[4]]
  
  boundaries <- (status_boundaries - RMA.int) / RMA.slope
  
  return(boundaries)
  
}

#' @describeIn boundaries Estimates nutrient boundaries using the "mismatch" method of minimising the difference between two types of missclassification.
#' @inheritParams boundaries_lm
#' @export
boundaries_MM <- function(data,
                          yvar = "y.u",
                          xvar = "x.u",
                          status_boundaries = c("HG" = 0.8, "GM" = 0.6, "MP" = 0.4)) {
  # Objective function
  f <- function(data, x, b) {
    Nutr.not.Good_EQR_Good = sum(data[[xvar]] > x &
                                   data[[yvar]] >= b, na.rm = T)
    Nutr.Good_EQR.not.Good = sum(data[[xvar]] < x &
                                   data[[yvar]] < b, na.rm = T)
    Mismatch_diff = (Nutr.not.Good_EQR_Good - Nutr.Good_EQR.not.Good) ^2
    return(Mismatch_diff)
  }
  
  boundaries <- lapply(status_boundaries, function(b)
    optimize(
      f,
      interval = c(min(data[[xvar]]), max(data[[xvar]])),
      data = data,
      b = b
    )$minimum)
  
  return(unlist(boundaries))
  
}

#' @describeIn boundaries Estimates nutrient boundaries using the "adjusted quartiles" method.
#' @inheritParams boundaries_lm
#' @param class_var name of variable defining status classes. Defaults to "status".
#' @export
boundaries_AdjQ <- function(data,
                          yvar = "y.u",
                          xvar = "x.u",
                          class_var = "status"
                          ) {

 q75s <- tapply(data[[xvar]], data[[class_var]], quantile, probs = 0.75)
 q25s <- tapply(data[[xvar]], data[[class_var]], quantile, probs = 0.25)
 
 HG <- mean(c(q75s[["High"]], q25s[["Good"]]))
 GM <- mean(c(q75s[["Good"]], q25s[["Moderate"]]))
 MP <- mean(c(q75s[["Moderate"]], q25s[["Poor"]]))
 
 return(c("HG"=HG, "GM"=GM, "MP"=MP))
  
}

#' @describeIn boundaries Estimates nutrient boundaries as the mean of the median concentrations in adjacent classes
#' @inheritParams boundaries_lm
#' @export
boundaries_medians <- function(data,
                            yvar = "y.u",
                            xvar = "x.u",
                            class_var = "status"
) {
  
  q5s <- tapply(data[[xvar]], data[[class_var]], quantile, probs = 0.5)
  
  HG <- mean(c(q5s[["High"]], q5s[["Good"]]))
  GM <- mean(c(q5s[["Good"]], q5s[["Moderate"]]))
  MP <- mean(c(q5s[["Moderate"]], q5s[["Poor"]]))
  
  return(c("HG"=HG, "GM"=GM, "MP"=MP))
  
}

#' @describeIn boundaries Estimates nutrient boundaries as the upper quartile of the relevant status class
#' @inheritParams boundaries_lm
#' @export
boundaries_q75 <- function(data,
                           yvar = "y.u",
                           xvar = "x.u",
                           class_var = "status"
) {
  
  q75s <- tapply(data[[xvar]], data[[class_var]], quantile, probs = 0.75)
  
  HG <- q75s[["High"]]
  GM <- q75s[["Good"]]
  MP <- q75s[["Moderate"]]
  
  return(c("HG"=HG, "GM"=GM, "MP"=MP))
  
}



# dat$status <- cut(dat$y.u, breaks = rev(c(Inf, status_boundaries, -Inf)), labels = rev(c("High", "Good", "Moderate", "Poor")))
# boundaries_q50(dat)
