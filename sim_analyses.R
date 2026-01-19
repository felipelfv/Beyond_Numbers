# COMPLETE ANALYSIS OF SIMULATION RESULTS AND PLOTS
library(ggplot2); library(tidyr); library(dplyr)

results <- simulation_results_factorial_with_nulls_100

# helper
parse_modelfit <- function(modelfit_text) {
  waic_line <- grep("WAIC", modelfit_text, value = TRUE)
  waic <- NA
  if (length(waic_line) > 0) {
    waic <- as.numeric(gsub(".*WAIC\\s+", "", waic_line[1]))
  }
  
  dic_line <- grep("DIC2", modelfit_text, value = TRUE)
  dic <- NA
  if (length(dic_line) > 0) {
    dic <- as.numeric(gsub(".*DIC2\\s+", "", dic_line[1]))
  }
  
  return(c(WAIC = waic, DIC = dic))
}

# ANALYSIS FUNCTION
analyze_condition <- function(condition_results, condition_name) {
  
  cat("CONDITION:", condition_name, "\n")
  
  n_reps <- length(condition_results)
  
  # DIF scenario
  dif_scenario <- condition_results[[1]]$dif_scenario
  dif_items <- condition_results[[1]]$dif_items
  
  cat("DIF items:", paste(dif_items, collapse = ", "), "\n")
  cat("Replications:", n_reps, "\n\n")
  
  # true values (note: same for both models)
  true_dif <- rep(0, 4)
  names(true_dif) <- c("y2", "y3", "y4", "y5")
  true_dif[dif_items] <- 0.30  # True DIF effect magnitude
  
  # DIF PARAMETER ESTIMATES - AGE MODEL
  age_dif_list <- list(y2 = list(), y3 = list(), y4 = list(), y5 = list())
  
  for (r in seq_along(condition_results)) {
    # to hndle both possible data structures
    if (!is.null(condition_results[[r]]$age$estimates)) {
      # old structure (without full/null)
      ests <- condition_results[[r]]$age$estimates
    } else {
      # new structure (with full/null)
      ests <- condition_results[[r]]$age$full$estimates
    }
    
    # standardized parameters
    age_dif_list$y2[[r]] <- ests["y2 ~ latX*lat", ]
    age_dif_list$y3[[r]] <- ests["y3 ~ latX*lat", ]
    age_dif_list$y4[[r]] <- ests["y4 ~ latX*lat", ]
    age_dif_list$y5[[r]] <- ests["y5 ~ latX*lat", ]
  }
  
  age_results <- list()
  
  for (item in c("y2", "y3", "y4", "y5")) {
    if (!is.null(age_dif_list[[item]][[1]])) {
      means <- sapply(age_dif_list[[item]], function(x) as.numeric(x["Mean"]))
      ci_lower <- sapply(age_dif_list[[item]], function(x) as.numeric(x["2.5%"]))
      ci_upper <- sapply(age_dif_list[[item]], function(x) as.numeric(x["97.5%"]))
      
      true_val <- true_dif[item]
      mean_est <- mean(means, na.rm = TRUE)
      bias <- mean_est - true_val
      rmse <- sqrt(mean((means - true_val)^2, na.rm = TRUE))
      coverage <- mean((true_val >= ci_lower) & (true_val <= ci_upper), na.rm = TRUE)
      
      age_results[[item]] <- list(
        mean = mean_est, 
        bias = bias, 
        rmse = rmse, 
        coverage = coverage
      )
    }
  }
  
  # DIF PARAMETER ESTIMATES - GROUP MODEL
  group_dif_list <- list(y2 = list(), y3 = list(), y4 = list(), y5 = list())
  
  for (r in seq_along(condition_results)) {
    if (!is.null(condition_results[[r]]$group$estimates)) {
      ests <- condition_results[[r]]$group$estimates
    } else {
      ests <- condition_results[[r]]$group$full$estimates
    }
    
    group_dif_list$y2[[r]] <- ests["y2 ~ Group*lat", ]
    group_dif_list$y3[[r]] <- ests["y3 ~ Group*lat", ]
    group_dif_list$y4[[r]] <- ests["y4 ~ Group*lat", ]
    group_dif_list$y5[[r]] <- ests["y5 ~ Group*lat", ]
  }
  
  group_results <- list()
  
  for (item in c("y2", "y3", "y4", "y5")) {
    if (!is.null(group_dif_list[[item]][[1]])) {
      means <- sapply(group_dif_list[[item]], function(x) as.numeric(x["Mean"]))
      
      # same true value as age model
      true_val <- true_dif[item]  # 0.30 for DIF items, 0 for non-DIF items
      mean_est <- mean(means, na.rm = TRUE)
      bias <- mean_est - true_val
      
      group_results[[item]] <- list(mean = mean_est, bias = bias)
    }
  }
  
  # WALD TEST REJECTION RATES
  # age model wald tests
  age_wald_pvals <- do.call(rbind, lapply(condition_results, function(x) {
    if (!is.null(x$age$wald)) {
      wald <- x$age$wald
    } else {
      wald <- x$age$full$wald
    }
    data.frame(
      b_int_y2 = wald[1, "probability"],
      b_int_y3 = wald[2, "probability"],
      b_int_y4 = wald[3, "probability"],
      b_int_y5 = wald[4, "probability"],
      joint = wald[5, "probability"]
    )
  }))
  
  age_power <- list()
  for (i in 1:4) {
    test_name <- c("b_int_y2", "b_int_y3", "b_int_y4", "b_int_y5")[i]
    reject_rate <- mean(age_wald_pvals[[test_name]] < 0.05, na.rm = TRUE)
    age_power[[test_name]] <- reject_rate
  }
  age_power$joint <- mean(age_wald_pvals$joint < 0.05, na.rm = TRUE)
  
  # group model wald tests
  group_wald_pvals <- do.call(rbind, lapply(condition_results, function(x) {
    if (!is.null(x$group$wald)) {
      wald <- x$group$wald
    } else {
      wald <- x$group$full$wald
    }
    data.frame(
      b_gint_y2 = wald[1, "probability"],
      b_gint_y3 = wald[2, "probability"],
      b_gint_y4 = wald[3, "probability"],
      b_gint_y5 = wald[4, "probability"],
      joint = wald[5, "probability"]
    )
  }))
  
  group_type1 <- list()
  for (test in colnames(group_wald_pvals)) {
    group_type1[[test]] <- mean(group_wald_pvals[[test]] < 0.05, na.rm = TRUE)
  }
  
  # MODEL FIT
  model_fit <- do.call(rbind, lapply(condition_results, function(x) {
    if (!is.null(x$group$modelfit)) {
      # old structure
      group_fit <- parse_modelfit(x$group$modelfit)
      age_fit <- parse_modelfit(x$age$modelfit)
      
      data.frame(
        group_waic = group_fit["WAIC"],
        group_dic = group_fit["DIC"],
        age_waic = age_fit["WAIC"],
        age_dic = age_fit["DIC"]
      )
    } else {
      # new structure with full/null
      group_full <- parse_modelfit(x$group$full$modelfit)
      group_null <- parse_modelfit(x$group$null$modelfit)
      age_full <- parse_modelfit(x$age$full$modelfit)
      age_null <- parse_modelfit(x$age$null$modelfit)
      
      data.frame(
        group_waic = group_full["WAIC"],
        group_dic = group_full["DIC"],
        age_waic = age_full["WAIC"],
        age_dic = age_full["DIC"],
        group_null_waic = group_null["WAIC"],
        group_null_dic = group_null["DIC"],
        age_null_waic = age_null["WAIC"],
        age_null_dic = age_null["DIC"]
      )
    }
  }))
  
  return(list(
    condition = condition_name,
    N = condition_results[[1]]$N,
    dif_scenario = dif_scenario,
    n_reps = n_reps,
    age_results = age_results,
    group_results = group_results,
    age_power = age_power,
    group_type1 = group_type1,
    model_fit = model_fit
  ))
}

all_summaries <- list()

for (condition_name in names(results)) {
  summary <- analyze_condition(results[[condition_name]], condition_name)
  all_summaries[[condition_name]] <- summary
}


# PLOTS
bias_colors <- c("eta_x" = "#2E7D32", "Group" = "#C62828")    
power_colors <- c("eta_x" = "#1565C0", "Group" = "#F57C00")   
fit_colors <- c("WAIC" = "#6A1B9A", "DIC" = "#00838F")      

bias_list <- list()
power_list <- list()
fit_list <- list()

for (cond in names(all_summaries)) {
  s <- all_summaries[[cond]]
  
  parts <- strsplit(cond, "_")[[1]]
  n_val <- as.numeric(gsub("N0*", "", parts[1]))
  dif_scenario <- parts[2]
  
  true_dif_items <- if (dif_scenario == "2DIF") c("y2", "y3") else c("y2", "y3", "y4", "y5")
  
  # eta_x model bias
  age_bias <- data.frame(
    Condition = cond,
    N = n_val,
    DIF_Scenario = dif_scenario,
    Model = "eta_x",
    Item = c("y2", "y3", "y4", "y5"),
    Bias = c(s$age_results$y2$bias, s$age_results$y3$bias, 
             s$age_results$y4$bias, s$age_results$y5$bias),
    True_DIF = c("y2", "y3", "y4", "y5") %in% true_dif_items
  )
  
  # group model bias - now using same true values
  group_bias <- data.frame(
    Condition = cond,
    N = n_val,
    DIF_Scenario = dif_scenario,
    Model = "Group",
    Item = c("y2", "y3", "y4", "y5"),
    Bias = c(s$group_results$y2$bias, s$group_results$y3$bias,
             s$group_results$y4$bias, s$group_results$y5$bias),
    True_DIF = c("y2", "y3", "y4", "y5") %in% true_dif_items
  )
  
  bias_list[[paste0(cond, "_age")]] <- age_bias
  bias_list[[paste0(cond, "_group")]] <- group_bias
  
  # power data
  power_list[[cond]] <- data.frame(
    Condition = cond,
    N = n_val,
    DIF_Scenario = dif_scenario,
    Item = c("y2", "y3", "y4", "y5"),
    eta_x_Power = c(s$age_power$b_int_y2, s$age_power$b_int_y3,
                    s$age_power$b_int_y4, s$age_power$b_int_y5),
    Group_Power = c(s$group_type1$b_gint_y2, s$group_type1$b_gint_y3,
                    s$group_type1$b_gint_y4, s$group_type1$b_gint_y5),
    True_DIF = c("y2", "y3", "y4", "y5") %in% true_dif_items
  )
  
  # model fit data
  fit_list[[cond]] <- data.frame(
    Condition = cond,
    N = n_val,
    DIF_Scenario = dif_scenario,
    WAIC_Diff = mean(s$model_fit$group_waic - s$model_fit$age_waic, na.rm = TRUE),
    DIC_Diff = mean(s$model_fit$group_dic - s$model_fit$age_dic, na.rm = TRUE)
  )
}

bias_data <- do.call(rbind, bias_list)
power_data <- do.call(rbind, power_list)
fit_data <- do.call(rbind, fit_list)

# BIAS (both models)
bias_data$N_label <- factor(paste0("N = ", bias_data$N), 
                            levels = paste0("N = ", c(200, 400, 800)))
bias_data$DIF_label <- factor(ifelse(bias_data$DIF_Scenario == "2DIF", 
                                     "2 DIF Items", "4 DIF Items"),
                              levels = c("2 DIF Items", "4 DIF Items"))

p1 <- ggplot(bias_data, aes(x = Item, y = Bias, fill = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_col(position = position_dodge(0.7), width = 0.7) +
  facet_grid(N_label ~ DIF_label) +
  scale_fill_manual(values = bias_colors, name = "Model") +
  labs(title = "", # true DIF = 0.30 for DIF items, 0 for non-DIF items
       x = "Item",
       y = "Bias") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"),
        panel.grid.minor = element_blank())

print(p1)
ggsave("plot1_dif_bias.pdf", p1, width = 8, height = 6)


# POWER/TYPE I ERROR
power_long <- power_data %>%
  pivot_longer(cols = c(eta_x_Power, Group_Power),
               names_to = "Model",
               values_to = "Reject_Rate") %>%
  mutate(Model = gsub("_Power", "", Model),
         N_label = factor(paste0("N = ", N), 
                          levels = paste0("N = ", c(200, 400, 800))),
         DIF_label = factor(ifelse(DIF_Scenario == "2DIF", 
                                   "2 DIF Items", "4 DIF Items"),
                            levels = c("2 DIF Items", "4 DIF Items")))

p2 <- ggplot(power_long, aes(x = Item, y = Reject_Rate, fill = Model)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_col(position = position_dodge(0.7), width = 0.7) +
  facet_grid(N_label ~ DIF_label) +
  scale_fill_manual(values = power_colors, name = "Model") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "", 
       subtitle = "",
       x = "Item",
       y = "Rejection Rate") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"),
        panel.grid.minor = element_blank())

print(p2)
ggsave("plot2_power_type1.pdf", p2, width = 8, height = 6)


# JOINT TEST POWER
joint_list <- list()
for (cond in names(all_summaries)) {
  s <- all_summaries[[cond]]
  parts <- strsplit(cond, "_")[[1]]
  n_val <- as.numeric(gsub("N0*", "", parts[1]))
  dif_scenario <- parts[2]
  
  joint_list[[cond]] <- data.frame(
    N = n_val,
    DIF_Scenario = dif_scenario,
    eta_x_Joint = s$age_power$joint,
    Group_Joint = s$group_type1$joint
  )
}

joint_data <- do.call(rbind, joint_list)

joint_long <- joint_data %>%
  pivot_longer(cols = c(eta_x_Joint, Group_Joint),
               names_to = "Model",
               values_to = "Reject_Rate") %>%
  mutate(Model = gsub("_Joint", "", Model),
         DIF_label = factor(ifelse(DIF_Scenario == "2DIF", 
                                   "2 DIF Items", "4 DIF Items"),
                            levels = c("2 DIF Items", "4 DIF Items")))

p3 <- ggplot(joint_long, aes(x = N, y = Reject_Rate, color = Model, shape = Model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  facet_wrap(~ DIF_label) +
  scale_color_manual(values = power_colors) +
  scale_shape_manual(values = c("eta_x" = 16, "Group" = 17)) +
  scale_x_continuous(breaks = c(200, 400, 800)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "", 
       x = "Sample Size",
       y = "Rejection Rate (Joint Wald Test)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"),
        panel.grid.minor = element_blank())

print(p3)
ggsave("plot3_joint_power.pdf", p3, width = 8, height = 6)

# MODEL FIT COMPARISON (WAIC & DIC)
fit_long <- fit_data %>%
  pivot_longer(cols = c(WAIC_Diff, DIC_Diff),
               names_to = "Criterion",
               values_to = "Difference") %>%
  mutate(
    Criterion = gsub("_Diff", "", Criterion),
    DIF_label = factor(ifelse(DIF_Scenario == "2DIF", 
                              "2 DIF Items", "4 DIF Items"),
                       levels = c("2 DIF Items", "4 DIF Items"))
  )

p4 <- ggplot(fit_long, aes(x = N, y = Difference, color = Criterion, shape = Criterion)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  facet_wrap(~ DIF_label) +
  scale_color_manual(values = fit_colors) +
  scale_shape_manual(values = c("WAIC" = 16, "DIC" = 17)) +
  scale_x_continuous(breaks = c(200, 400, 800)) +
  labs(title = "",  
       subtitle = "", 
       x = "Sample Size",
       y = "Mean Difference in Criterion") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"),
        panel.grid.minor = element_blank())

print(p4)
ggsave("plot4_model_fit.pdf", p4, width = 8, height = 6)

# FULL VS NULL MODEL
if ("group_null_waic" %in% colnames(all_summaries[[1]]$model_fit)) {
  
  full_null_list <- list()
  
  for (cond in names(all_summaries)) {
    s <- all_summaries[[cond]]
    parts <- strsplit(cond, "_")[[1]]
    n_val <- as.numeric(gsub("N0*", "", parts[1]))
    dif_scenario <- parts[2]
    
    full_null_list[[cond]] <- data.frame(
      N = n_val,
      DIF_Scenario = dif_scenario,
      Group_WAIC = mean(s$model_fit$group_null_waic - s$model_fit$group_waic, na.rm = TRUE),  
      Group_DIC = mean(s$model_fit$group_null_dic - s$model_fit$group_dic, na.rm = TRUE),      
      eta_x_WAIC = mean(s$model_fit$age_null_waic - s$model_fit$age_waic, na.rm = TRUE),     
      eta_x_DIC = mean(s$model_fit$age_null_dic - s$model_fit$age_dic, na.rm = TRUE)        
    )
  }
  
  full_null_data <- do.call(rbind, full_null_list)
  
  full_null_long <- full_null_data %>%
    pivot_longer(cols = c(Group_WAIC, Group_DIC, eta_x_WAIC, eta_x_DIC),
                 names_to = "Metric",
                 values_to = "Diff") %>%
    mutate(
      Model = if_else(grepl("Group", Metric), "Group", "eta_x"),
      Criterion = if_else(grepl("WAIC", Metric), "WAIC", "DIC"),
      DIF_label = factor(ifelse(DIF_Scenario == "2DIF", 
                                "2 DIF Items", "4 DIF Items"),
                         levels = c("2 DIF Items", "4 DIF Items"))
    )
  
  p5 <- ggplot(full_null_long, aes(x = N, y = Diff, color = Model, linetype = Criterion)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~ DIF_label) +
    scale_color_manual(values = c("eta_x" = unname(fit_colors["WAIC"]), 
                                  "Group" = unname(fit_colors["DIC"]))) +
    scale_linetype_manual(values = c("WAIC" = "solid", "DIC" = "twodash")) +
    scale_x_continuous(breaks = c(200, 400, 800)) +
    labs(title = "", 
         subtitle = "", 
         x = "Sample Size",
         y = "Mean Difference in Criterion") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  print(p5)
  ggsave("plot5_full_vs_null.pdf", p5, width = 8, height = 6)
}

saveRDS(all_summaries, "analysis_summaries.rds")