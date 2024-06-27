# library("fixest")
# library("plm")

#' Pre-Mean Differencing
#'
#' Calculate pre-treatment means rather than just using the y value from t-1.
#' @param df The dataset used in the analysis.
#' @param y The character string that denotes the outcome variable.
#' @param pmd_lag The number of pre-periods that should be used when calculating the pre-treatment mean.
#' @return The reference value in the analysis.
#' @noRd
pmd_func <- function(df, y, pmd_lag){

  # This function averages y-values from t-1 to t-pmd_lag
  the_lag <- 0
  for(k in 1:pmd_lag) the_lag <- the_lag + lag(df[,y], k)
  the_lag <- the_lag / pmd_lag
  return(the_lag)
}

#' Post-Mean Differencing
#'
#' Calculate post-treatment means rather than just using a single y value from t+h.
#' @param df The dataset used in the analysis.
#' @param y The character string that denotes the outcome variable.
#' @param pool_lead The number of post-periods that should be used when calculating the post-treatment mean.
#' @return The average of all future values in the analysis.
#' @noRd
pooled_adjustment <- function(df, y, pool_lead){

  the_lead <- 0
  for(k in 0:pool_lead) the_lead <- the_lead + lead(df[,y], k)
  the_lead <- the_lead / (pool_lead+1)
  return(the_lead)
}

#' Get LP-DiD Weights
#'
#' Calculate weights to use in the post-treatment periods.
#' @param df The dataset used in the analysis.
#' @param j The variable indicating the current horizon.
#' @param time_index The variable indicating calendar time.
#' @return Regression weights for the post-treatment horizon periods.
#' @noRd
get_weights <- function(df, j, time_index, lim){

  df[,paste0("group_h", j)] <- NA
  df[,paste0("group_h", j)] <- ifelse(lim, df[,time_index], df[,paste0("group_h", j)])

  frmla <- as.formula(paste0("treat_diff ~ 1 | ", time_index))
  tmp <- feols(frmla, data = df[lim,])

  df[,paste0("num_weights_", j)] <- NA
  df[lim,paste0("num_weights_", j)] <- tmp$residuals
  df[is.na(df$treat_diff) | df$treat_diff == 0,paste0("num_weights_", j)] <- NA
  den_weights = sum(df[,paste0("num_weights_", j)], na.rm = T)
  df[,paste0("weight_", j)] <- df[,paste0("num_weights_", j)] / den_weights
  suppressWarnings(df[,paste0("gweight_", j)] <- ave(df[,paste0("weight_", j)], df[,paste0("group_h", j)], FUN = function(x) max(x, na.rm = TRUE)))

  lim <- is.na(df[,paste0("weight_", j)])
  df[lim,paste0("weight_", j)] <- df[lim,paste0("gweight_", j)]

  tmp <- 1 / df[,paste0("weight_", j)]
  df <- df[,colnames(df)[!grepl("group_|weights_|weight_[^0]|.weight_", colnames(df))]]
  return(tmp)
}

#' Plot LP-DiD Event Study
#'
#' Plotting LP-DiD event study parameter estimates and confidence intervals.
#' @param reg An object generated via the lpdid function.
#' @param conf The confidence level (1-alpha) desired for confidence intervals. Default is 0.95 (95% confidence internval)
#' @param segments A boolean (TRUE or FALSE) value for whether confidence intervals should be generated.  Default is TRUE.
#' @param add A boolean (TRUE or FALSE) value for whether the user wants to add to an existing figure or generate a new one.  The default is FALSE, meaning a new plot.
#' @param xlab The text belonging on the x axis.
#' @param ylab The text belonging on the y axis.
#' @param ylim A vector of length two that determines the minimum and maximum of the y-axis. Default scales to fit the figure being plotted.
#' @param main The text belonging as the figure's title.
#' @param x.shift A numeric value that will shift the event study estimates along the x-axis.
#' @param pch A numeric value corresponding to the point's shape.
#' @param cex A numeric value corresponding to the size of the point.
#' @param col The color the points (and confidence intervals) should be.
#' @param opacity A numeric value between 0 and 1 that corresponds to the opacity of the color.
#' @return An event study plot.
#' @export
plot_lpdid <- function(reg, conf = .95, segments = TRUE, add = FALSE,
                       xlab = NULL, ylab = NULL,
                       ylim = NULL,
                       main = "", x.shift = 0,
                       pch = 19, cex = 1, col = "black", opacity = 1){

  if(nrow(reg$coeftable) != length(reg$window)) stop("coeftable and window are not the same length.  It is likely that pooled=TRUE in the lpdid function.  An event study cannot be plotted when pooled=TRUE.")
  coeftable <- reg$coeftable
  coeftable$t <- reg$window
  conf_z <- abs(qnorm((1-conf)/2))
  uCI <- coeftable$Estimate + conf_z*coeftable$`Std. Error`
  lCI <- coeftable$Estimate - conf_z*coeftable$`Std. Error`

  if(!add){

    if(is.null(ylim)){

      ylim <- c(min(lCI, na.rm = T), max(uCI, na.rm = T))
    }

    plot(coeftable$t + x.shift, coeftable$Estimate, las = 1, pch = pch,
         cex = cex, main = main, col = scales::alpha(col, opacity),
         ylim = ylim,
         xlim = range(coeftable$t, na.rm = T),
         xlab = ifelse(is.null(xlab), "Time to Treatment", xlab),
         ylab = ifelse(is.null(ylab),
                       paste0("Coefficient Estimate and ", conf*100, "% Confidence Interval"),
                       ylab))
    abline(h = 0, v = -1, lty = 2)
  } else {

    points(coeftable$t + x.shift, coeftable$Estimate, pch = pch,
           cex = cex, col = scales::alpha(col, opacity))
  }

  if(segments) segments(x0 = coeftable$t + x.shift,
                        y0 = uCI, y1 = lCI, col = scales::alpha(col, opacity))
}

#' Simulate Data from Dube et al. (2023)
#'
#' This function generates simulated data from Dube et al. (2023).
#' @param exogenous_timing A boolean (TRUE or FALSE) value that, if set to TRUE, will make treatment exogenous.  When set to FALSE,
#'  previous values of Y will make treatment more likely, or endogenous. Default is TRUE
#' @param seed A numeric value to set.seed().
#' @return A simulated dataset.
#' @export
genr_sim_data <- function(exogenous_timing = TRUE, seed = 757){

  ## Set parameters ####
  N <- 500; T <- 50; e_sd <- 25; rho <- 0.5; alpha0 <- 2; alpha1 <- 0.5
  psi <- 0.6
  K <- 3
  # exogenous_timing <- T
  post_window <- 10
  pre_window <- 5

  set.seed(seed)

  ## Create Base Data ####

  df <- data.frame(i = rep(1:N, each = T),
                   t = rep(1:T, N),
                   e = rnorm(N*T, 0, e_sd))

  df <- pdata.frame(df, index=c("i","t"), drop.index=FALSE, row.names=FALSE)
  df$i <- as.numeric(df$i)
  df$t <- as.numeric(df$t)

  mu_unit <- rnorm(N, 0, e_sd); mu_year <- rnorm(T, 0, e_sd)
  df$y <- mu_unit[df$i] + mu_year[df$t] + df$e
  df$y <- df$y + ifelse(is.na(lag(df$y, 1)), 0, lag(df$y, 1)) # y_{it}(0)

  ## Exogenous Treatment ####

  if(exogenous_timing){

    groupz <- sample(rep(1:10, each = N/10), N, replace = FALSE)
    treat_datez <- seq(11, 27, by = 2)

    df$treat_date <- treat_datez[groupz[df$i]]
    df$treat_status <- ifelse(is.na(df$treat_date), 0, ifelse(df$t == df$treat_date, 1, 0))
    df$ever_treat <- ifelse(!is.na(df$treat_date), 1, 0)
    df$rel_time <- ifelse(df$ever_treat == 1, df$t - df$treat_date, df$t - (max(df$t) + 1))

    beta <- ifelse(df$rel_time < 0, 0,
                   ifelse(df$rel_time <= 20,
                          (alpha0 * (df$rel_time+1)) + (alpha1 * (df$rel_time+1)^2) + ((1-alpha1)*((df$rel_time+1)^2 / (df$treat_date / min(treat_datez))^2)),
                          (alpha0 * 21) + (alpha1 * 21^2) + ((1-alpha1)*(21^2 / (df$treat_date / min(treat_datez))^2))))
    df$y <- df$y + beta
  }

  ## Endogenous Treatment ####

  if(!exogenous_timing){

    df$y_diff <- df$y - lag(df$y, 1)
    df$sigma <- ave(df$y_diff, df$t, FUN = sd)

    cond1 <- (psi * lag(df$y_diff, 1)) + ((1-psi)*(rnorm(nrow(df), 0, e_sd))) < -1 * df$sigma
    cond2 <- df$t %in% 11:30
    g <- aggregate(df$t, list(ifelse(cond1 & cond2, 1, 0), df$i), min, drop = F)
    g <- g[g$Group.1 == 1,]
    table(g$x, useNA = "always")

    treat_datez <- sort(unique(g$x))

    df$treat_date <- g$x[match(df$i, g$Group.2)]
    df$treat_status <- ifelse(is.na(df$treat_date), 0, ifelse(df$t == df$treat_date, 1, 0))
    df$ever_treat <- ifelse(!is.na(df$treat_date), 1, 0)
    df$rel_time <- ifelse(df$ever_treat == 1, df$t - df$treat_date, df$t - (max(df$t) + 1))

    beta <- ifelse(df$rel_time < 0, 0,
                   ifelse(df$rel_time <= 20,
                          (alpha0 * (df$rel_time+1)) + (alpha1 * (df$rel_time+1)^2) + ((1-alpha1)*((df$rel_time+1)^2 / (df$treat_date / min(treat_datez))^2)),
                          (alpha0 * 21) + (alpha1 * 21^2) + ((1-alpha1)*(21^2 / (df$treat_date / min(treat_datez))^2))))
    df$y <- df$y + beta
  }

  df$beta <- beta
  return(df)
}

#' Local Projections Difference-in-Differences
#'
#' This function estimates LP-DiD regressions as outlined in Dube et al. (2023) <doi:10.3386/w31184>.
#' @param df The dataset used in the analysis.
#' @param window A vector of length two that denotes the length of the pre- and post-periods.
#'  The first number, denoting the number of pre-periods before treamtent, should be negative and second, denoting the number of post periods after treatment, should be positive.
#' @param depvar The outcome variable.  This should be input as a character, the name of the outcome variable in the data. Note: @param y will continue to work for backwards compatibility.
#' @param treat_status The name of the column that defines treatment. This column is expected to contain 0's before treatment, and 1's during and following treatment.
#' @param unit The name of the column that represents the unit ID.  This should be a character. Note: @param unit_index will continue to work for backwards compatibility.
#' @param time The name of the column that represents the calendar time.  This should be a character. Note: @param time_index will continue to work for backwards compatibility.
#' @param cluster The name of the column by which to cluster the standard errors.  Default is the unit ID.
#' @param controls A formula of control variables to be included in the regression formula. The form should be: `~ X1 + X2 | FE1 + FE2`. These controls will be time invariant.
#' @param controls_t A formula of control variables to be included in the regression formula. The form should be: `~ X1 + X2 | FE1 + FE2`. These controls will be time varying.
#' @param dylags The number of outcome lags to be included in the analysis.
#'  For an example of this, simulate endogenous data via genr_sim_data(FALSE), and include one outcome lag. Note: @param outcome_lags will continue to work for backwards compatibility.
#' @param rw A boolean (TRUE or FALSE) value that will re-weight the regression and generate ATT rather than VWATT. The default is FALSE, which corresponds to the estimator calculating the VWATT (variance weighted average treatment effect on the treated). Note: @param reweight will continue to work for backwards compatibility.
#' @param pmd A boolean (TRUE or FALSE) value that, if equal to TRUE, will use pre-treatment means rather than a single value from t-1.
#' @param pmd_lag The number of pre-treatment periods to include in taking the mean.
#' @param nocomp A boolean value that will remove later-treated observations from the control group even before they are treated.  See Section 2.10 "Composition effects". Note: @param composition_correction will continue to work for backwards compatibility.
#' @param pooled A boolean value (TRUE or FALSE) that, if equal to TRUE, will calculate the treatment effect pooled over all post-periods.
#' @param nonabsorbing_lag Sets the number of periods after which dynamic effects are assumed to stabilize (and thus the unit can return to the control group). This number must be greater than zero, and will indicate to the function that you'd like to estimate a nonabsorbing treatment effect.
#' @return A list including a coefficient table and window vector.
#' @export
lpdid <- function(df, window = c(NA, NA),
                  depvar,
                  unit, time,
                  treat_status = "",
                  cluster = NULL,
                  controls = NULL,
                  controls_t = NULL,
                  dylags = 0,
                  rw = FALSE,
                  pmd = FALSE, pmd_lag,
                  nocomp = FALSE,
                  pooled = FALSE,
                  nonabsorbing_lag = NULL,
                  y, unit_index, time_index,
                  outcome_lags = 0, reweight = FALSE,
                  composition_correction = FALSE){

  if(is.null(cluster)) cluster <- unit_index

  y <- ifelse(!missing(depvar), depvar, y)
  unit_index <- ifelse(!missing(unit), unit, unit_index)
  time_index <- ifelse(!missing(time), time, time_index)
  outcome_lags <- ifelse(!missing(dylags), dylags, outcome_lags)
  reweight <- ifelse(!missing(rw), rw, reweight)
  composition_correction <- ifelse(!missing(nocomp), nocomp, composition_correction)

  pre_window <- -1*window[1]; post_window <- window[2]

  if(!inherits(df, "pdata.frame")) df <- pdata.frame(df, index=c(unit_index,time_index), drop.index=FALSE, row.names=FALSE)
  df[,unit_index] <- as.character(df[,unit_index]); df[,time_index] <- as.numeric(df[,time_index])

  # Calculate Pre-Mean Differences
  if(pmd){

    df$the_lag <- pmd_func(df = df, y = y, pmd_lag = pmd_lag)
  } else {

    df$the_lag <- lag(df[,y], 1)
  }

  nonabsorbing <- F
  if(is.null(nonabsorbing_lag)) nonabsorbing_lag <- Inf else nonabsorbing <- T

  # create "rel_time" variable for nonabsorbing
  rel_time <- "rel_time"
  df[,rel_time] <- NA; df[,"treat_date"] <- NA

  # @param treat_status The name of the column that denotes treatment status.
  #  This vector should take on a value of 1 for each time the unit is treated.
  #  As an example, in a state-year panel, if a state is treated in 1990 and 2010, the state should have 1's for 1990 and 2010 only.  All other year observations of this variable for this unit should be equal to 0.

  if(FALSE){ # this is a placeholder for more complicated treatments...
    for(i in unique(df[,unit_index])){

      id_logic <- which(df[,unit_index] == i)
      valz <- df[id_logic[df[id_logic, treat_status] == 1], time_index]
      if(!is.na(valz[1])){

        mat <- matrix(NA, ncol = length(valz), nrow = length(id_logic))

        for(j in valz){

          mat[,which(j == valz)] <- df[id_logic, time_index] - j
        }

        df[id_logic, rel_time] <- apply(mat, 1, function(x) x[which(abs(x) == min(abs(x)))[1]])
        df[id_logic, "treat_date"] <- apply(mat, 1, function(x) valz[which(abs(x) == min(abs(x)))[1]])
      } else {

        df[id_logic, rel_time] <- -1000
      }
    }
  } else {

    agg <- aggregate(list(time_index = df[df[,treat_status] > 0,time_index]),
                     list(unit_index = df[df[,treat_status] > 0,unit_index]),
                     min)
    m <- match(df[,unit_index], agg$unit_index)
    df[,"treat_date"] <- ifelse(is.na(m), -1000, agg$time_index[m])
    df[,rel_time] <- ifelse(is.na(m), -1000, df[,time_index] - agg$time_index[m])
  }

  df$treat <- ifelse(df[,rel_time] >= 0, 1, 0)
  df$treat_diff <- df$treat - lag(df$treat, 1)
  df$treat_diff <- ifelse(df$treat_diff < 0, 0, df$treat_diff)

  # ## Insert Controls here...
  # controls <- ~ a + b | DV3
  # ## Insert Time-Varying Controls here...
  # controls_t <- ~ x + y + z
  # controls_t <- ~ x + y + z | FE1 + FE2

  controls <- as.character(controls)[2]
  controls_t <- as.character(controls_t)[2]

  controls <- strsplit(controls, "\\s*\\|\\s*")[[1]]
  if(length(controls) > 2 | length(controls_t) > 2) stop("An unexpected number of vertical bars included in one of the controls arguments.")

  FE <- NULL
  if(length(controls) == 2){

    FE <- controls[2]
    controls <- controls[1]
  }

  controls_t <- strsplit(controls_t, "\\s*\\|\\s*")[[1]]
  if(length(controls_t) == 2){

    if(is.null(FE)) FE <- controls_t[2] else FE <- c(FE, controls_t[2])
    controls_t <- controls_t[1]
  }

  controls <- strsplit(controls, "\\s*\\+\\s*")[[1]]
  controls_t <- strsplit(controls_t, "\\s*\\+\\s*")[[1]]

  # Calculate lags of the outcome
  lagz <- c()
  if(outcome_lags>0){

    for(outcome_lag in 1:outcome_lags){

      df[,paste0("y_diff_lag", outcome_lag)] <- lag(df[,y] - lag(df[,y], 1), outcome_lag)
    }
    lagz <- c(lagz, colnames(df)[grepl("y_diff_lag.", colnames(df))])
  }

  lpdid_betaz <- rep(0, length(-pre_window:post_window))
  lpdid_sez <- rep(0, length(-pre_window:post_window))
  lpdid_nz <- rep(0, length(-pre_window:post_window))
  lpdid_regz <- vector(mode = "list", length(-pre_window:post_window))

  if(pooled) loop_bound <- 0 else loop_bound <- max(post_window, pre_window)
  # j<-2
  for(j in 0:loop_bound){

    # Calculate Pooled
    if(pooled & j == 0) df[,y] <- pooled_adjustment(df, y, post_window)

    if(!reweight) df$reweight_0 <- df$reweight_use <- 1

    # Post
    if(j <= post_window){

      df$Dy <- lead(df[,y], j) - df$the_lag

      frmla <- "Dy ~ treat_diff"
      if(!is.null(lagz[1])) lagz <- paste(lagz, collapse = " + ")
      if(!is.null(lagz[1])) frmla <- paste0(frmla, " + ", lagz)
      if(!is.na(controls[1])) controls_use <- paste(controls, collapse = " + ")
      if(!is.na(controls[1])) frmla <- paste(frmla, "+", paste(controls, collapse = " + "))

      if(!is.na(controls_t[1])){

        for(ctrl in controls_t){

          df[,paste0(ctrl, ".l")] <- lead(df[,ctrl], j)
        }

        controls_t_use <- paste0(controls_t, ".l")
        frmla <- paste0(frmla, " + ", paste(controls_t_use, collapse = " + "))
      }

      if(!is.null(FE)){
        rhs <- paste(time_index, " + ", paste(unlist(strsplit(FE, "\\s*\\+\\s*")), collapse = " + "))
      } else {
        rhs <- time_index
      }
      frmla <- as.formula(paste0(frmla, " | ", rhs))

      df$cluster_var <- df[,cluster]

      # Create "Limit"
      if(!nonabsorbing){

        lim <- !is.na(df[,"Dy"]) & !is.na(df$treat_diff) & !is.na(lead(df$treat, j))
        if(composition_correction){

          lim <- lim & (df$treat_diff == 1 | lead(df$treat, post_window) == 0) & (is.na(df$treat_date) | (df$treat_date < max(df[,time_index]) - post_window))
        } else {

          lim <- lim & (df$treat_diff == 1 | lead(df$treat, j) == 0)
        }
      } else {

        ## Non-Absorbing Limit
        lim_ctrl <- TRUE; lim_treat <- TRUE
        for(i in -nonabsorbing_lag:j){

          lim_ctrl <- lim_ctrl & lag(df$treat_diff, -i) == 0
          lim_treat <- lim_treat & if(i >= 0) lead(df$treat, i) == 1 else lag(df$treat, -i) == 0
        }
        df$lim_c <- ifelse(lim_ctrl, 1, 0)
        df$lim_t <- ifelse(lim_treat, 1, 0)
        lim <- lim_ctrl | lim_treat
      }
      lim <- !is.na(lim) & lim

      # Calculate weights for j
      if(reweight){

        df$reweight_use <- get_weights(df = df, j = j, time_index = time_index, lim = lim)
        if(j == 0) df$reweight_0 <- df$reweight_use
      }

      lim <- lim & df$reweight_use > 0

      # Check for perfect multi-colinearity
      if(sum(lim) > 0 && sum(!aggregate(df$treat_diff[lim], list(df[lim, time_index]), mean)$x %in% c(0, 1))){

        # Estimate and Save
        tmp <- suppressMessages(feols(frmla, data = df[lim,], cluster = ~cluster_var, weights = ~reweight_use))
        lpdid_betaz[match(j, -pre_window:post_window)] <- tmp$coeftable[1,1]
        lpdid_sez[match(j, -pre_window:post_window)] <- tmp$coeftable[1,2]
        lpdid_nz[match(j, -pre_window:post_window)] <- nobs(tmp)
        lpdid_regz[[match(j, -pre_window:post_window)]] <- tmp
      } else {

        lpdid_betaz[match(j, -pre_window:post_window)] <- NA
        lpdid_sez[match(j, -pre_window:post_window)] <- NA
        lpdid_nz[match(j, -pre_window:post_window)] <- NA
      }
    }

    # Pre
    if((j>1 & j<=pre_window) || (pmd & j<=pre_window)){

      df$Dy <- lag(df[,y], j) - df$the_lag

      frmla <- "Dy ~ treat_diff"
      if(!is.null(lagz[1])) lagz <- paste(lagz, collapse = " + ")
      if(!is.null(lagz[1])) frmla <- paste0(frmla, " + ", lagz)
      if(!is.na(controls[1])) controls_use <- paste(controls, collapse = " + ")
      if(!is.na(controls[1])) frmla <- paste(frmla, "+", paste(controls, collapse = " + "))

      if(!is.na(controls_t[1])){

        for(ctrl in controls_t){

          df[,paste0(ctrl, ".l")] <- lag(df[,ctrl], j)
        }

        controls_t_use <- paste0(controls_t, ".l")
        frmla <- paste0(frmla, " + ", paste(controls_t_use, collapse = " + "))
      }

      if(!is.null(FE)){
        rhs <- paste(time_index, " + ", paste(unlist(strsplit(FE, "\\s*\\+\\s*")), collapse = " + "))
      } else {
        rhs <- time_index
      }
      frmla <- as.formula(paste0(frmla, " | ", rhs))

      df$cluster_var <- df[,cluster]

      if(!nonabsorbing){

        lim <- !is.na(df[,"Dy"]) & !is.na(df$treat_diff) & !is.na(df$treat)

        if(composition_correction){
          lim <- lim & (df$treat_diff == 1 | lead(df$treat, post_window) == 0) & (is.na(df$treat_date) | (df$treat_date < (max(df[,time_index]) - post_window)))
        } else {
          lim <- lim & (df$treat_diff == 1 | df$treat == 0)
        }

      } else {

        ## Non-Absorbing Limit
        lim_ctrl <- TRUE; lim_treat <- TRUE
        for(i in -nonabsorbing_lag:j){

          lim_ctrl <- lim_ctrl & lag(df$treat_diff, -i) == 0
          lim_treat <- lim_treat & if(i >= 0) lead(df$treat, i) == 1 else lag(df$treat, -i) == 0
        }
        lim <- lim_ctrl | lim_treat
      }

      lim <- !is.na(lim) & lim & df$reweight_0 > 0

      if(sum(lim) > 0 && sum(!aggregate(df$treat_diff[lim], list(df[lim, time_index]), mean)$x %in% c(0, 1))){

        suppressMessages(tmp <- feols(frmla, data = df[lim,], cluster = ~cluster_var, weights = ~reweight_0))
        lpdid_betaz[match(-j, -pre_window:post_window)] <- tmp$coeftable[1,1]
        lpdid_sez[match(-j, -pre_window:post_window)] <- tmp$coeftable[1,2]
        lpdid_nz[match(-j, -pre_window:post_window)] <- nobs(tmp)
        lpdid_regz[[match(-j, -pre_window:post_window)]] <- tmp
      } else {

        lpdid_betaz[match(-j, -pre_window:post_window)] <- NA
        lpdid_sez[match(-j, -pre_window:post_window)] <- NA
        lpdid_nz[match(-j, -pre_window:post_window)] <- NA
      }
    }
  }

  coeftable <- data.frame(Estimate = lpdid_betaz,
                          "Std. Error" = lpdid_sez,
                          "t value" = lpdid_betaz/lpdid_sez,
                          "Pr(>|t|)" = NA,
                          check.names = FALSE)
  if(pooled) coeftable <- coeftable[match(0, -pre_window:post_window),]
  coeftable[,4] <- pnorm(abs(coeftable$`t value`), lower.tail = F)

  print(coeftable[!is.na(coeftable$Estimate),])
  invisible(list(coeftable = coeftable[!is.na(coeftable$Estimate),],
                 # df = df,
                 window = c(-pre_window:post_window)[!is.na(coeftable$Estimate)],
                 nobs = data.frame(window = c(-pre_window:post_window),
                                   nobs = lpdid_nz),
                 regz = lpdid_regz))
}
