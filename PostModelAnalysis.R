# R version
# R 4.5.2

# Attached packages and versions
# - broom.mixed     (v0.2.9.7)
# - emmeans         (v2.0.0)
# - GGally          (v2.4.0)
# - ggplot2         (v4.0.1)
# - glmmTMB         (v1.1.15.2)
# - latex2exp       (v0.9.6)
# - lme4            (v1.1.37)
# - mgcv            (v1.9.3)
# - multcompView    (v0.1.10)
# - patchwork       (v1.3.2)
# - tidyverse       (v2.0.0)
# - vioplot         (v0.5.1)


library(lme4)         # Linear and generalized linear mixed-effects models
library(emmeans)      # Estimated marginal means and post hoc comparisons
library(tidyverse)    # Data manipulation, visualization, and wrangling tools
library(mgcv)         # Generalized additive models (GAMs)
library(glmmTMB)      # Flexible GLMMs (Beta, Gamma, Poisson, etc.)
library(purrr)        # Functional programming and iteration tools
library(ggplot2)      # Data visualization and plotting
library(latex2exp)    # Render LaTeX mathematical expressions in ggplot
library(broom.mixed)  # Convert mixed-model outputs into tidy data frames
library(vioplot)      # Violin plots for visualizing distributions
library(patchwork)    # Combine multiple ggplots into a single figure
library(multcompView) # Compact letter displays for multiple comparisons
library(GGally)       # Extensions to ggplot2 (e.g., pairwise plots)



# Load the 500 simulated networks
sim_networks <- readRDS("~/Trait-mediated_mutualistic_network/sim_networks.rds")
data <- sim_networks



# ==========================================================
#### Effects of adaptations on guild-level outcomes and network structure ####
# ==========================================================
#
# Scenarios:
#   Baseline : ecological dynamics only
#   Beh      : behavioural adaptation
#   Evo      : trait evolution
#   Evo-Beh  : evolution and behaviour combined
#
# For each response variable, a GLMM is fitted and treatment
# contrasts are calculated relative to the ecological-only
# scenario ("Eco").
#
# ==========================================================

# ----------------------------------------------------------
# Function: Calculate treatment contrasts
# ----------------------------------------------------------

contrast_fun <- function(var, data) {

  dat <- data %>%
    select(!any_of(ends_with("0"))) %>%
    select(starts_with(var)) %>%
    set_names(c(
      "Eco+Evo+Beh",
      "Eco+Evo",
      "Eco+Beh",
      "Eco"
    )) %>%
    mutate(rep = as.integer(rownames(.))) %>%
    pivot_longer(
      cols = c(
        "Eco+Evo+Beh",
        "Eco+Evo",
        "Eco+Beh",
        "Eco"
      ),
      names_to = "scenario",
      values_to = "response"
    ) %>%
    filter(is.finite(response))

  # Beta and Gamma distributions require positive values
  if (grepl("rao|modularity", var)) {
    dat <- dat %>%
      filter(response > 0)
  }

  # --------------------------------------------------------
  # Select model family
  # --------------------------------------------------------

  fam <- case_when(
    grepl("connectance|modularity", var) ~ "beta",
    grepl("persis|sum|rao", var) ~ "gamma",
    grepl("nestedness", var) ~ "gaussian"
  )

  # --------------------------------------------------------
  # Fit GLMM
  # --------------------------------------------------------

  m <- switch(
    fam,

    beta = glmmTMB(
      response ~ scenario + (1 | rep),
      data = dat,
      family = beta_family(link = "logit")
    ),

    gamma = glmmTMB(
      response ~ scenario + (1 | rep),
      data = dat,
      family = Gamma(link = "log")
    ),

    gaussian = glmmTMB(
      response ~ scenario + (1 | rep),
      data = dat,
      family = gaussian(link = "identity")
    )
  )

  # --------------------------------------------------------
  # Estimated marginal means and contrasts
  # --------------------------------------------------------

  emm <- emmeans(m, ~ scenario)

  contrasts <- as.data.frame(
    contrast(
      emm,
      method = "trt.vs.ctrl",
      ref = "Eco"
    )
  )

  contrasts$contrast <- c("Beh", "Evo", "Evo-Beh")
  contrasts$var <- var

  contrasts
}


# ----------------------------------------------------------
# Function: Build contrast data frame
# ----------------------------------------------------------

contrast_df_fun <- function(var_set, data) {

  contrast_df <- list_rbind(
    purrr::map(var_set, ~ contrast_fun(.x, data))
  )

  contrast_df <- contrast_df %>%
    mutate(

      response = case_when(
        grepl("persis", var) ~ "persistence",
        grepl("sum", var) ~ "biomass",
        grepl("rao", var) ~ "functional diversity (Rao's)",
        TRUE ~ as.character(var)
      ),

      guild = case_when(
        stringr::str_sub(var, -2) %in% c("_P", "_x") ~ "plant",
        stringr::str_sub(var, -2) %in% c("_A", "_y") ~ "pollinator",
        TRUE ~ "network"
      ),

      sig = case_when(
        p.value <= 0.05 ~ "*",
        p.value > 0.05 ~ "",
        TRUE ~ NA_character_
      )
    )

  contrast_df$response <- factor(
    contrast_df$response,
    levels = c(
      "persistence",
      "functional diversity (Rao's)",
      "biomass",
      "weighted connectance",
      "modularity Q",
      "weighted nestedness"
    )
  )

  contrast_df$guild <- factor(
    contrast_df$guild,
    levels = c(
      "plant",
      "pollinator",
      "network"
    )
  )

  contrast_df
}


# ----------------------------------------------------------
# Function: Plot treatment effects
# ----------------------------------------------------------

plot_contrast <- function(df) {

  ggplot(df, aes(
    x = contrast,
    y = estimate,
    fill = guild
  )) +

    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +

    geom_col(
      position = position_dodge(width = 0.7),
      alpha = 0.8,
      width = 0.6
    ) +

    geom_errorbar(
      aes(
        ymin = estimate - 1.96 * SE,
        ymax = estimate + 1.96 * SE
      ),
      width = 0.2,
      linewidth = 0.8,
      colour = "black",
      position = position_dodge(width = 0.7)
    ) +

    geom_text(
      aes(
        label = sig,
        y = ifelse(
          estimate >= 0,
          estimate + 1.96 * SE,
          estimate - 1.96 * SE
        )
      ),
      position = position_dodge(width = 0.7),
      vjust = ifelse(df$estimate >= 0, -0.1, 1.3),
      size = 6,
      colour = "black"
    ) +

    scale_fill_manual(
      values = c(
        plant = "forestgreen",
        pollinator = "goldenrod",
        network = "skyblue"
      )
    ) +

    facet_wrap(
      ~ response,
      ncol = 2,
      scales = "free_y"
    ) +

    scale_y_continuous(
      limits = function(x) {
        c(
          x[1] + x[1] * 0.10,
          x[2] + x[2] * 0.10
        )
      }
    ) +

    labs(
      x = "Scenario",
      y = "Effect of adaptation",
      fill = "Guild"
    ) +

    theme_minimal(base_size = 16)
}


# ----------------------------------------------------------
# Variables analysed
# ----------------------------------------------------------

var_set <- c(
  "persis_P",
  "persis_A",
  "sum_P",
  "sum_A",
  "rao_P",
  "rao_A",
  "weighted connectance",
  "modularity Q",
  "weighted nestedness"
)

# ----------------------------------------------------------
# Run analysis
# ----------------------------------------------------------

contrast_df <- contrast_df_fun(
  var_set = var_set,
  data = data
)

contrast_plot <- plot_contrast(contrast_df)

# Display figure
contrast_plot

# ==========================================================
#### FIGURE 4 Effects of Kernel Widths Across Adaptation Scenarios ####
# ==========================================================
#
# This analysis evaluates how variation in kernel widths
# (σm, σk, σc) influences community and network responses
# under four adaptation scenarios:

#
# Fixed effects represent changes in response variables with
# increasing kernel width. Interaction terms quantify
# scenario-specific deviations from the ecological baseline.
#
# ==========================================================

# ----------------------------------------------------------
# Response variable labels
# ----------------------------------------------------------

response_labels <- c(
  persis_P = "(a) Plant persistence",
  persis_A = "(b) Pollinator persistence",
  sum_P    = "(c) Plant biomass",
  sum_A    = "(d) Pollinator biomass",
  rao_P    = "(e) Plant functional diversity",
  rao_A    = "(f) Pollinator functional diversity",
  `weighted connectance` = "(g) Weighted connectance",
  `modularity Q` = "(h) Modularity",
  `weighted nestedness` = "(i) Nestedness"
)

# ----------------------------------------------------------
# Extract kernel-width effects for a response variable
# ----------------------------------------------------------

kernel_fun <- function(var, data) {

  dat <- data %>%
    select(!any_of(ends_with("0"))) %>%
    select(starts_with(var), c("km", "kk", "kc")) %>%
    set_names(c(
      "Both",
      "Evo",
      "Beh",
      "baseline",
      "km",
      "kk",
      "kc"
    )) %>%
    mutate(rep = as.integer(rownames(.))) %>%
    pivot_longer(
      cols = c("Both", "Evo", "Beh", "baseline"),
      names_to = "scenario",
      values_to = "response"
    ) %>%
    filter(is.finite(response)) %>%
    mutate(
      exp_km = exp(km),
      exp_kk = exp(kk),
      exp_kc = exp(kc),
      scenario = factor(scenario),
      rep = factor(rep)
    )

  # Positive-valued distributions require positive responses
  if (grepl("rao|modularity", var)) {
    dat <- dat %>%
      filter(response > 0)
  }

  # --------------------------------------------------------
  # Select model family
  # --------------------------------------------------------

  fam <- case_when(
    grepl("connectance|modularity", var) ~ "beta",
    grepl("persis|sum|rao", var) ~ "gamma",
    grepl("nestedness", var) ~ "gaussian"
  )

  # --------------------------------------------------------
  # Fit GLMM
  # --------------------------------------------------------

  m <- switch(

    fam,

    beta = glmmTMB(
      response ~
        (exp_km + exp_kk + exp_kc) * scenario +
        (1 | rep),
      family = beta_family(link = "logit"),
      data = dat
    ),

    gamma = glmmTMB(
      response ~
        (exp_km + exp_kk + exp_kc) * scenario +
        (1 | rep),
      family = Gamma(link = "log"),
      data = dat
    ),

    gaussian = glmmTMB(
      response ~
        (exp_km + exp_kk + exp_kc) * scenario +
        (1 | rep),
      family = gaussian(link = "identity"),
      data = dat
    )
  )

  # --------------------------------------------------------
  # Extract coefficients
  # --------------------------------------------------------

  coef_df <- tidy(m, effects = "fixed") %>%
    select(!any_of(c("effect", "component"))) %>%
    filter(grepl("_", term))

  # Baseline slopes
  km_slope <- coef_df$estimate[1]
  kk_slope <- coef_df$estimate[2]
  kc_slope <- coef_df$estimate[3]

  # --------------------------------------------------------
  # Compute scenario-specific slopes
  # --------------------------------------------------------

  coef_df <- coef_df %>%
    mutate(

      total_estimate = case_when(
        !grepl(":", term) ~ estimate,
        grepl("km", term) ~ estimate + km_slope,
        grepl("kk", term) ~ estimate + kk_slope,
        grepl("kc", term) ~ estimate + kc_slope
      ),

      conf.low  = total_estimate - 1.96 * std.error,
      conf.high = total_estimate + 1.96 * std.error,

      var = var,

      kernel = case_when(
        grepl("km", term) ~ "km",
        grepl("kk", term) ~ "kk",
        grepl("kc", term) ~ "kc",
        TRUE ~ NA_character_
      ),

      scenario = case_when(
        grepl("Beh", term) ~ "Beh",
        grepl("Evo", term) ~ "Evo",
        grepl("Both", term) ~ "Evo-Beh",
        TRUE ~ "None"
      ),

      term = case_when(
        grepl("km", term) ~ "exp_km",
        grepl("kk", term) ~ "exp_kk",
        grepl("kc", term) ~ "exp_kc",
        TRUE ~ NA_character_
      )
    )

  coef_df
}

# ----------------------------------------------------------
# Run analysis across all response variables
# ----------------------------------------------------------

kernel_df <- list_rbind(
  purrr::map(
    var_set,
    ~ kernel_fun(.x, data)
  )
)

# ----------------------------------------------------------
# Format plotting data
# ----------------------------------------------------------

kernel_df <- kernel_df %>%
  mutate(

    var = factor(
      var,
      levels = names(response_labels),
      labels = response_labels
    ),

    scenario = factor(
      scenario,
      levels = c(
        "None",
        "Beh",
        "Evo",
        "Evo-Beh"
      )
    ),

    sig = if_else(p.value <= 0.05, "*", "")
  )

# ----------------------------------------------------------
# Plot kernel-width effects
# ----------------------------------------------------------

kernel_plot <- ggplot(
  kernel_df,
  aes(
    x = total_estimate,
    y = term,
    xmin = conf.low,
    xmax = conf.high,
    colour = scenario,
    shape = scenario
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey50"
  ) +

  geom_errorbar(
    position = position_dodge(width = 0.6),
    width = 0.2
  ) +

  geom_point(
    position = position_dodge(width = 0.6),
    size = 2.5
  ) +

  scale_shape_manual(
    values = c(
      "None"    = 1,
      "Beh"     = 17,
      "Evo"     = 15,
      "Evo-Beh" = 16
    )
  ) +

  scale_y_discrete(
    labels = c(
      "exp_kk" = TeX("$\\sigma_k$"),
      "exp_kc" = TeX("$\\sigma_c$"),
      "exp_km" = TeX("$\\sigma_m$")
    )
  ) +

  facet_wrap(
    ~ var,
    ncol = 2,
    scales = "free_x"
  ) +

  labs(
    x = "Estimate (95% CI)",
    y = "Kernel width",
    colour = "Scenario",
    shape = "Scenario"
  ) +

  theme_classic(base_size = 16)

# Display figure
kernel_plot

# ==========================================================
#### FIGURE 5. Effects of Adaptive Foraging and Trait Evolution on Trait-Mediated Coefficients ####
# ==========================================================
#
# This analysis compares the distributions of trait-mediated
# coefficients across four adaptation scenarios:
#
#   Baseline : ecological dynamics only
#   Beh      : behavioural adaptation
#   Evo      : trait evolution
#   Evo-Beh  : evolution and behaviour combined
#
# For each coefficient, linear mixed-effects models are used
# to test for differences among scenarios. Estimated marginal
# means are compared using Sidak-adjusted post hoc tests and
# displayed with compact letter displays.
#
# ==========================================================

# ----------------------------------------------------------
# Function: Violin plot with mixed-model comparisons
# ----------------------------------------------------------

plot_violin <- function(
    data,
    vars,
    labels,
    fill_col,
    title,
    ylab = "Mean"
) {

  df <- data %>%
    mutate(rep = factor(seq_len(nrow(.)))) %>%
    select(all_of(c("rep", vars))) %>%
    pivot_longer(
      cols = all_of(vars),
      names_to = "Scenario",
      values_to = "Value"
    ) %>%
    mutate(
      Scenario = factor(
        Scenario,
        levels = vars,
        labels = labels
      )
    )

  # --------------------------------------------------------
  # Linear mixed-effects model
  # --------------------------------------------------------

  mod <- lmer(
    Value ~ Scenario + (1 | rep),
    data = df
  )

  # --------------------------------------------------------
  # Estimated marginal means
  # --------------------------------------------------------

  emm <- emmeans(mod, ~ Scenario)

  # --------------------------------------------------------
  # Compact letter display
  # --------------------------------------------------------

  cld_df <- as.data.frame(
    multcomp::cld(
      emm,
      adjust = "sidak",
      Letters = letters
    )
  )

  letter_df <- data.frame(
    Scenario = cld_df$Scenario,
    Letters = cld_df$.group,
    stringsAsFactors = FALSE
  )

  # Position letters above violins

  ypos <- df %>%
    group_by(Scenario) %>%
    summarise(
      y = max(Value) * 1.30,
      .groups = "drop"
    )

  letter_df <- left_join(
    letter_df,
    ypos,
    by = "Scenario"
  )

  # --------------------------------------------------------
  # Plot
  # --------------------------------------------------------

  ggplot(
    df,
    aes(
      x = Scenario,
      y = Value
    )
  ) +

    geom_violin(
      fill = fill_col,
      alpha = 0.8,
      trim = FALSE
    ) +

    geom_boxplot(
      width = 0.10,
      fill = "white",
      outlier.shape = NA
    ) +

    geom_text(
      data = letter_df,
      aes(
        x = Scenario,
        y = y,
        label = Letters
      ),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    ) +

    labs(
      title = title,
      x = NULL,
      y = ylab
    ) +

    theme_classic(base_size = 14) +

    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

# ----------------------------------------------------------
# Scenario labels
# ----------------------------------------------------------

labels <- c(
  "Baseline",
  "Beh",
  "Evo",
  "Evo-Beh"
)

# ==========================================================
# Mean coefficient distributions
# ==========================================================

p1 <- plot_violin(
  data = data,
  vars = c(
    "mean_cPf3",
    "mean_cPf2",
    "mean_cPf1",
    "mean_cPf"
  ),
  labels = labels,
  fill_col = "forestgreen",
  title = "Plant\ncompetition coefficient"
)

p2 <- plot_violin(
  data = data,
  vars = c(
    "mean_KPf3",
    "mean_KPf2",
    "mean_KPf1",
    "mean_KPf"
  ),
  labels = labels,
  fill_col = "forestgreen",
  title = "Plant\ncarrying capacity"
)

p3 <- plot_violin(
  data = data,
  vars = c(
    "mean_cAf3",
    "mean_cAf2",
    "mean_cAf1",
    "mean_cAf"
  ),
  labels = labels,
  fill_col = "goldenrod",
  title = "Pollinator\ncompetition coefficient"
)

p4 <- plot_violin(
  data = data,
  vars = c(
    "mean_KAf3",
    "mean_KAf2",
    "mean_KAf1",
    "mean_KAf"
  ),
  labels = labels,
  fill_col = "goldenrod",
  title = "Pollinator\ncarrying capacity"
)

p5 <- plot_violin(
  data = data,
  vars = c(
    "mean_sigmaf3",
    "mean_sigmaf2",
    "mean_sigmaf1",
    "mean_sigmaf"
  ),
  labels = labels,
  fill_col = "skyblue",
  title = "Plant/Pollinator\nbenefit coefficient"
)

# ==========================================================
# Standard deviation distributions
# ==========================================================

p6 <- plot_violin(
  data = data,
  vars = c(
    "sd_cPf3",
    "sd_cPf2",
    "sd_cPf1",
    "sd_cPf"
  ),
  labels = labels,
  fill_col = "forestgreen",
  title = "Plant\ncompetition coefficient",
  ylab = "SD"
)

p7 <- plot_violin(
  data = data,
  vars = c(
    "sd_KPf3",
    "sd_KPf2",
    "sd_KPf1",
    "sd_KPf"
  ),
  labels = labels,
  fill_col = "forestgreen",
  title = "Plant\ncarrying capacity",
  ylab = "SD"
)

p8 <- plot_violin(
  data = data,
  vars = c(
    "sd_cAf3",
    "sd_cAf2",
    "sd_cAf1",
    "sd_cAf"
  ),
  labels = labels,
  fill_col = "goldenrod",
  title = "Pollinator\ncompetition coefficient",
  ylab = "SD"
)

p9 <- plot_violin(
  data = data,
  vars = c(
    "sd_KAf3",
    "sd_KAf2",
    "sd_KAf1",
    "sd_KAf"
  ),
  labels = labels,
  fill_col = "goldenrod",
  title = "Pollinator\ncarrying capacity",
  ylab = "SD"
)

p10 <- plot_violin(
  data = data,
  vars = c(
    "sd_sigmaf3",
    "sd_sigmaf2",
    "sd_sigmaf1",
    "sd_sigmaf"
  ),
  labels = labels,
  fill_col = "skyblue",
  title = "Plant/Pollinator\nbenefit coefficient",
  ylab = "SD"
)

# ----------------------------------------------------------
# Combine panels
# ----------------------------------------------------------

violin_plot <-
  (p1 | p2 | p3 | p4 | p5) /
  (p6 | p7 | p8 | p9 | p10)

# Display figure

violin_plot

# ----------------------------------------------------------
# Save figure
# ----------------------------------------------------------

ggsave(
  filename = "violin_plot.jpeg",
  plot = violin_plot,
  width = 13,
  height = 7,
  dpi = 400
)

